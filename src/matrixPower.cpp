#include "matrixPower.h"
#include "utility.h"
#include <cmath>
#include "macros.h"

mtxPower::mtxPower(Graph* graph_, int highestPower_, double cacheSize_, double safetyFactor_):graph(graph_),levelPtr(NULL),macroLevelPtr(NULL),levelData(NULL), highestPower(highestPower_), cacheSize(cacheSize_), safetyFactor(safetyFactor_)
{
    traverser = new Traverse(graph, RACE::ONE);
}

mtxPower::~mtxPower()
{
    if(levelPtr)
    {
        delete[] levelPtr;
    }
    if(macroLevelPtr)
    {
        delete[] macroLevelPtr;
    }
    delete traverser;
    if(levelData)
    {
        delete levelData;
    }
}

double mtxPower::getElemUpperLimit(int level)
{
    return (safetyFactor*(2*highestPower+1)*levelData->levelNnz[level]);
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
    printf("cacheElem = %f, nnzr = %f\n", cacheElem, nnzr);

    std::vector<int> cacheViolatedLevel;
    std::vector<int> cacheViolatedFactor;
    //check if there is a level where nnz violates cache
    //so the one where it violates first is detected first
    for(int level=0; level<totalLevel; ++level)
    {
        //printf("NNZ[%d] = %d\n", level, levelData->levelNnz[level]);
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
    consolidatePartition();
}

void mtxPower::consolidatePartition()
{
    std::vector<int> newLevelPtr;
    newLevelPtr.push_back(0);
    double sumElem = 0;
    int sumNnz = 0;
    int ctr = 0;
    for(int level=0; level<totalLevel; ++level)
    {
        double currElem = getElemUpperLimit(level);
        double bytePerNNZ = getBytePerNNZ();
        double cacheElem = cacheSize/bytePerNNZ;
        sumElem += currElem;
        sumNnz += levelData->levelNnz[level];
        if(sumElem >= cacheElem)
        {
            sumElem = currElem;
            newLevelPtr.push_back(levelPtr[level]);
            //rewrite levelRow and levelNnz
            int curr_idx = (int)newLevelPtr.size();
            levelData->levelRow[ctr] = newLevelPtr[curr_idx-1] - newLevelPtr[curr_idx-2];
            levelData->levelNnz[ctr] = (sumNnz-levelData->levelNnz[level]);
            sumNnz = levelData->levelNnz[level];
            ++ctr;
        }
    }
    newLevelPtr.push_back(levelPtr[totalLevel]);
    int curr_idx = (int)newLevelPtr.size();
    levelData->levelRow[ctr] = newLevelPtr[curr_idx-1] - newLevelPtr[curr_idx-2];
    levelData->levelNnz[ctr] = sumNnz;

    //rewrite levelPtr with newLevelPtr
    totalLevel = (int)(newLevelPtr.size()-1);
    printf("consolidated levels = %d\n", totalLevel);
    for(int i=0; i<totalLevel+1; ++i)
    {
        levelPtr[i] = newLevelPtr[i];
    }

    for(int level=0; level<totalLevel; ++level)
    {
        printf("ROW[%d] = %d, NNZ[%d] = %d\n", level, levelData->levelRow[level], level, levelData->levelNnz[level]);
    }

#if 0 //for debugging
    printf("Validation of levelNnz and levelRow\n");
    printf("totalLevel = %d, ctr = %d\n", totalLevel, ctr);

    int sumRow = 0;
    sumNnz = 0;
    for(int i=0; i<totalLevel; ++i)
    {
        sumRow += levelData->levelRow[i];
        sumNnz += levelData->levelNnz[i];
    }
    printf("Total rows = %d, total nnz = %d\n", sumRow, sumNnz);
#endif

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
