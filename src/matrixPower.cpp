#include "matrixPower.h"
#include "utility.h"
#include <cmath>
#include "macros.h"
#include "lb.h"

mtxPower::mtxPower(Graph* graph_, int highestPower_, int numSharedCache_, double cacheSize_, double safetyFactor_):graph(graph_),levelPtr(NULL),levelGroupPtr(NULL),cacheLevelGroup(NULL),levelData(NULL), highestPower(highestPower_), numSharedCache(numSharedCache_), cacheSize(cacheSize_), safetyFactor(safetyFactor_)
{
    traverser = new Traverse(graph, RACE::ONE);
}

mtxPower::~mtxPower()
{
    if(levelPtr)
    {
        delete[] levelPtr;
    }
    if(levelGroupPtr)
    {
        delete[] levelGroupPtr;
    }
    if(cacheLevelGroup)
    {
        delete[] cacheLevelGroup;
    }
    delete traverser;
    if(levelData)
    {
        delete levelData;
    }
}

double mtxPower::getElemUpperLimit(int level)
{
    return (safetyFactor*(2*highestPower-1)*levelData->levelNnz[level]);
}

void mtxPower::createLevelPtr()
{
    levelPtr = new int[totalLevel+1];

    levelPtr[0] = 0;

    for(int level=0; level<totalLevel; ++level)
    {
        levelPtr[level+1] = levelPtr[level] + levelData->levelRow[level];
    }
}

double mtxPower::getBytePerNNZ()
{
    double nnzr = levelData->nnz/(double)levelData->nrow;
    //(1+power/nnzr) -> includes mtx values, lhs and rhs vector values
    //(1+1/nnzr) -> includes col and rowPtr
    return (double)(((1+highestPower/nnzr)*sizeof(double)+(1+1/nnzr)*sizeof(int)));
}

void mtxPower::findPartition()
{
    traverser->calculateDistance();
    levelData = traverser->getLevelData();
    totalLevel = levelData->totalLevel;
    createLevelPtr();
    double nnzr = levelData->nnz/(double)levelData->nrow;
    //convert cache Size to elements,
    //i.e. nnz for this matrix that cache can hold
    double bytePerNNZ = getBytePerNNZ();
    double cacheElem = cacheSize/bytePerNNZ;
    printf("cacheElem = %f\n", cacheElem);
    printf("nrows = %d\n", levelData->nrow);
    printf("nnz = %d\n", levelData->nnz);
    printf("nnzr = %f\n", nnzr);

    std::vector<int> cacheViolatedLevel;
    std::vector<int> cacheViolatedFactor;
    //check if there is a level where nnz violates cache
    //so the one where it violates first is detected first
    for(int level=0; level<totalLevel; ++level)
    {
        printf("NNZ[%d] = %d\n", level, levelData->levelNnz[level]);
        double currElem = levelData->levelNnz[level];
        //depending on power we have to adapt cache condition
        /*for(int idx=0; idx<highestPower; ++idx)
        {
            currElem += levelData->levelNnz[level-idx];
        }*/

        currElem = getElemUpperLimit(level);
        //if violated mark the levels
        if(currElem > cacheElem)
        {
            cacheViolatedLevel.push_back(level);
            int factor = static_cast<int>(ceil(getElemUpperLimit(level)/cacheElem));
            cacheViolatedFactor.push_back(factor);
            printf("Cache Violated at %d level, elem = %f, factor = %d\n", level, currElem, factor);
        }
    }

#if 0 //recursion step
    //try splitting by (factor+1) parts the first level that violates cache
    //(factor+1) parts since it while again doing BFS it may incur parts in
    //same level, therefore not able to maintain exactly the factor
    long int rowsToSplit = levelPtr[cacheViolatedLevel[0]+1] - levelPtr[cacheViolatedLevel[0]];
    long int subRowSize = static_cast<long int> (rowsToSplit/(double)cacheViolatedFactor[0]);
    int range_lo = levelPtr[cacheViolatedLevel[0]];
    int range_hi = graph->NROW; //we check from range_lo to entire range
    printf("lo = %d, hi = %d\n", range_lo, range_hi);
    //now traverse within first part
    Traverse* subTraverser = new Traverse(graph, RACE::ONE, range_lo, range_hi, 0, subRowSize);
    subTraverser->calculateDistance();
    LevelData* subLevelData = subTraverser->getLevelData();

    for(int level=0; level<subLevelData->totalLevel; ++level)
    {
        printf("NNZ[%d] = %d\n", level, subLevelData->levelNnz[level]);
    }
#endif
    splitSharedCacheDomain();
    consolidatePartition();
}

//split into 'n' equal parts, where 'n' is number of shared caches available
void mtxPower::splitSharedCacheDomain()
{
    std::vector<double> eff_vec;
    getEnv("RACE_EFFICIENCY", eff_vec);

    if(eff_vec.empty())
    {
        eff_vec.push_back(50);
    }
    //use lb for load balancing for nSharedCache
    LB lb(numSharedCache, eff_vec[0], levelData, RACE::POWER);
    lb.balance();
    int len;
    lb.getZonePtr(&cacheLevelGroup, &len);

    if(numSharedCache != (len-1))
    {
        WARNING_PRINT("All cache groups cannot be active, active  = %d, requested = %d, check reducing efficiency", len-1, numSharedCache);
        numSharedCache = (len-1);
    }
    /*printf("caches = %d, len = %d\n", numSharedCache, len);
    for(int i=0; i<(len-1); ++i)
    {
        printf("lg[%d] = %d, lg[%d] = %d, row = %d\n", i, cacheLevelGroup[i], i+1, cacheLevelGroup[i+1], cacheLevelGroup[i+1]-cacheLevelGroup[i]);
    }*/

    levelGroupPtr = new int[len];
    findMacroLevelPtr(cacheLevelGroup, levelGroupPtr);

    for(int i=0; i<len; ++i)
    {
        printf("lg = %d\n", levelGroupPtr[i]);
    }
}

void mtxPower::findMacroLevelPtr(int* zones, int* macroLevelPtr)
{
    int ctr=0;
    for(int i=0; i<(totalLevel+1); ++i)
    {
        if(levelPtr[i] == zones[ctr])
        {
            macroLevelPtr[ctr] = i;
            ++ctr;
        }
    }
}

void mtxPower::consolidatePartition()
{
    std::vector<int> newLevelPtr;
    newLevelPtr.push_back(0);
    for(int levelGroup=0; levelGroup<numSharedCache; ++levelGroup)
    {
        double sumElem = 0;
        for(int level=levelGroupPtr[levelGroup]; level<levelGroupPtr[levelGroup+1]; ++level)
        {
            double currElem = getElemUpperLimit(level);
            double bytePerNNZ = getBytePerNNZ();
            double cacheElem = cacheSize/bytePerNNZ;
            sumElem += currElem;
            if(sumElem >= cacheElem)
            {
                sumElem = currElem;
                newLevelPtr.push_back(levelPtr[level]);
            }
        }
        int lastLevel = levelGroupPtr[levelGroup+1];
        newLevelPtr.push_back(levelPtr[lastLevel]);
    }

    //rewrite levelPtr with newLevelPtr
    totalLevel = (int)(newLevelPtr.size()-1);
    printf("consolidated levels = %d\n", totalLevel);
    for(int i=0; i<totalLevel+1; ++i)
    {
        levelPtr[i] = newLevelPtr[i];
        printf("levelPtr[%d] = %d\n", i, levelPtr[i]);
    }

    int* newLevelGroupPtr = new int[numSharedCache+1];
    findMacroLevelPtr(cacheLevelGroup, newLevelGroupPtr);

    for(int i=0; i<(numSharedCache+1); ++i)
    {
        levelGroupPtr[i] = newLevelGroupPtr[i];
        printf("lg_consolidated[%d] = %d, row = %d\n", i, levelGroupPtr[i], levelPtr[levelGroupPtr[i]]);
    }

    delete[] newLevelGroupPtr;
}

void mtxPower::getPerm(int **perm_, int *len)
{
    traverser->getPerm(perm_, len);
}

void mtxPower::getInvPerm(int **invPerm_, int *len)
{
    traverser->getInvPerm(invPerm_, len);
}

int mtxPower::getTotalLevel()
{
    return totalLevel;
}

LevelData* mtxPower::getLevelDataRef()
{
    return levelData;
}

int* mtxPower::getLevelPtrRef()
{
    return levelPtr;
}

int mtxPower::getTotalLevelGroup()
{
    return numSharedCache;
}

int* mtxPower::getLevelGroupPtrRef()
{
    return levelGroupPtr;
}
