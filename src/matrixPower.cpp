#include "matrixPower.h"
#include "utility.h"
#include <cmath>
#include "macros.h"
#include "lb.h"
#include "omp.h"

//#define LB_REMINDER

mtxPower::mtxPower(Graph* graph_, int highestPower_, int numSharedCache_, double cacheSize_, double safetyFactor_):graph(graph_), levelPtr(NULL),levelGroupPtr(NULL), unlockRow(NULL), dangerRow(NULL), unlockCtr(NULL), cacheLevelGroup(NULL),levelData(NULL), highestPower(highestPower_), numSharedCache(numSharedCache_), cacheSize(cacheSize_), safetyFactor(safetyFactor_)
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
    if(unlockRow)
    {
        delete[] unlockRow;
    }
    if(dangerRow)
    {
        delete[] dangerRow;
    }
    if(unlockCtr)
    {
        delete[] unlockCtr;
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
    //(highestPower+1)*NNZ
    return (safetyFactor*(highestPower+1)*levelData->levelNnz[level]);
    //return (safetyFactor*(2*highestPower-1)*levelData->levelNnz[level]);
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

std::vector<int> identifyHopelessRegions(std::vector<int> cacheViolatedLevel)
{
    if(cacheViolatedLevel.size() > 1)
    {
        cacheViolatedLevel.push_back(cacheViolatedLevel.back()+2); //a dummy so last element will be pushed, the dummy is made 2 apart from last so it doesn't compress
    }
    std::vector<int> hopelessRegion;
    std::vector<int> hopelessRegion_flag(cacheViolatedLevel.size(), 0);
    int prevHopelessLevel = -2;
    int hopelessCtr = 0;
    int hopelessBlockCtr = 1;
    for(int i=0; i<(int)cacheViolatedLevel.size(); ++i)
    {
        int currHopelessLevel = cacheViolatedLevel[i];
        //printf("currHopeless = %d\n", currHopelessLevel);
        if(currHopelessLevel == (prevHopelessLevel+1))
        {
            //if consequtive
            hopelessRegion_flag[i-1] = hopelessBlockCtr;
            hopelessRegion_flag[i] = hopelessBlockCtr;
            hopelessCtr++;
        }
        else
        {
            hopelessBlockCtr++;
        }

        prevHopelessLevel = currHopelessLevel;
    }

   /* for(int i=0; i<(int)cacheViolatedLevel.size(); ++i)
    {
        printf("hopelessRegion_flag[%d] = %d\n", i, hopelessRegion_flag[i]);
    }*/
    printf("Compressing %d hopeless levels\n", hopelessCtr);
    int prevFlag = 0;
    //now compress parts that have the same flag set
    for(int i=0; i<(int)hopelessRegion_flag.size(); ++i)
    {
        int currFlag = hopelessRegion_flag[i];
        if(currFlag != prevFlag)
        {
            if((currFlag > 0) && (prevFlag == 0))
            {
                //this is start of a hopeless region
                hopelessRegion.push_back(cacheViolatedLevel[i]);
            }
            else if((currFlag > 0) && (prevFlag > 0))
            {
                hopelessRegion.push_back(cacheViolatedLevel[i-1]);
                hopelessRegion.push_back(cacheViolatedLevel[i]);
            }
            else if(currFlag == 0)
            {
                //this is end of a hopeless region
                hopelessRegion.push_back(cacheViolatedLevel[i-1]);
            }

            //printf("@@@@ check hopelessRegion[%d] = %d\n", (int)hopelessRegion.size()-1, hopelessRegion.back());
        }
        prevFlag = currFlag;
    }
    return hopelessRegion;
}

int mtxPower::get_cache_violation_cutoff(int stage)
{
    if(stage < (int)cache_violation_cutoff.size())
    {
        return cache_violation_cutoff[stage];
    }
    else
    {
        return cache_violation_cutoff[0];
    }
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

    int default_cutoff = 10; //make 10 factor so not much cut-off happens
    cache_violation_cutoff.push_back(default_cutoff); //default
    getEnv("RACE_CACHE_VIOLATION_CUTOFF", cache_violation_cutoff);
    printf("RACE_CACHE_VIOLATION_CUTOFF = %d\n", get_cache_violation_cutoff(1));
    for(int level=0; level<totalLevel; ++level)
    {
        printf("rowStart = %d, NROW[%d] = %d, NNZ[%d] = %d\n", levelPtr[level], level, levelData->levelRow[level], level, levelData->levelNnz[level]);
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
            int factor = static_cast<int>(ceil(getElemUpperLimit(level)/cacheElem));

            //for merging, merge nearby levels that have violation greater than
            //cacheViolationTolerance, because they are not going to bring
            //anything, then why spending expensive barrier cost
            if(factor > get_cache_violation_cutoff(1))
            {
                cacheViolatedLevel.push_back(level);
                cacheViolatedFactor.push_back(factor);
            }
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
    std::vector<int> hopelessRegion = identifyHopelessRegions(cacheViolatedLevel);
    consolidatePartition(hopelessRegion);
    findUnlockCtr();
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

    getStatNUMA();
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

void getHopelessStartEnd(std::vector<int> hopelessRegion, int count, int *start, int *end)
{
    (*start) = -1;
    (*end) = -1;
    if(2*count <  (int)hopelessRegion.size()-1)
    {
        (*start) = hopelessRegion[2*count];
        (*end) = hopelessRegion[2*count+1];
    }
}

void mtxPower::consolidatePartition(std::vector<int> hopelessRegion)
{
    unlockRow = new int[totalLevel];
    dangerRow = new int[totalLevel];

    std::vector<int> newLevelPtr;
    newLevelPtr.push_back(0); //done in loop
    int consolidated_ctr = 0;
    unlockRow[0] = levelPtr[1];

#ifdef LB_REMINDER
    //find max start and end boundary levels
    std::vector<int> maxSumBoundaries(highestPower-1,0); //counter for each boundary level

    for(int pow=0; pow<highestPower-1; ++pow)
    {
        for(int levelGroup=0; levelGroup<numSharedCache; ++levelGroup)
        {
            int startBoundary = levelGroupPtr[levelGroup]+pow;
            int endBoundary = levelGroupPtr[levelGroup+1]-1-pow;
            int startNNZ = levelData->levelNnz[startBoundary];
            int endNNZ =  levelData->levelNnz[endBoundary];

            maxSumBoundaries[pow] = std::max(maxSumBoundaries[pow],startNNZ+endNNZ);
        }
        printf("Boundary @ power %d = %d\n", pow+1, maxSumBoundaries[pow]);
    }
#endif

    int hopelessRegionCount = 0;
    int currHopelessStart, currHopelessEnd;
    getHopelessStartEnd(hopelessRegion, hopelessRegionCount, &currHopelessStart, &currHopelessEnd);

    for(int levelGroup=0; levelGroup<numSharedCache; ++levelGroup)
    {
        double sumElem = 0;
        double sumNNZ = 0;
        int curLevelCtr = 0;
        int consolidated_curLevelCtr = 0;
        //int totalLevelInGroup = levelGroupPtr[levelGroup+1]-levelGroupPtr[levelGroup];


        //newLevelPtr.push_back(levelPtr[levelGroupPtr[levelGroup]]);
        for(int level=levelGroupPtr[levelGroup]; level<levelGroupPtr[levelGroup+1]; ++level)
        {
            double currElem = getElemUpperLimit(level);
            double bytePerNNZ = getBytePerNNZ();
            double cacheElem = cacheSize/bytePerNNZ;
            sumElem += currElem;

            double curNNZ = levelData->levelNnz[level];
            sumNNZ += curNNZ;

            bool hopeless = false;
            //deal with hopeless region
            if((level >= currHopelessStart) && (level < currHopelessEnd))
            {
                //this level belongs to hopeless region
                hopeless = true;
            }
            else if(level == currHopelessEnd)
            {
                //don't set hopeless because now we have to end
                //and consolidate
                hopelessRegionCount++;
                getHopelessStartEnd(hopelessRegion, hopelessRegionCount, &currHopelessStart, &currHopelessEnd);
            }

            //printf("##### hopeless = %d, hopeless = [%d,%d]\n", hopeless, currHopelessStart, currHopelessEnd);
#ifdef LB_REMINDER
            int startPow = 0;
#endif
            if(curLevelCtr > 0 && !hopeless)
            {
                //second condition to avoid bulk level group at boundary,
                //as this is processed without caching

#ifdef LB_REMINDER
                if( (sumElem >= cacheElem) || ( (numSharedCache > 1) && ( (consolidated_curLevelCtr < (highestPower)) || (curLevelCtr > (totalLevelInGroup-highestPower)) ) ) )
#else
                if(sumElem >= cacheElem)
#endif
                {
                    bool updateFlag = true;

#ifdef LB_REMINDER
                    int endBoundary = levelGroupPtr[levelGroup+1]-(startPow)-1;
                    int endNNZ = levelData->levelNnz[endBoundary];

                    //for start boundary check if it has satisfied sum condition,
                    //not doing for end boundary since it is tricky there, as we are
                    //starting from start boundary and we would need to unroll back.
                    //But since only sum contributes to the load before barrier in
                    //reminder loop it should be fine.
                    if(consolidated_curLevelCtr < (highestPower))
                    {
                        updateFlag = false;
                        if((sumNNZ+endNNZ) > 0.85*maxSumBoundaries[startPow]) //0.85 for safety
                        {
                            updateFlag = true;
                            ++startPow;
                        }
                    }
#endif

                    if(updateFlag)
                    {
                        sumElem = currElem;
                        sumNNZ = curNNZ;
                        newLevelPtr.push_back(levelPtr[level]);

                        dangerRow[consolidated_ctr] = levelPtr[level-1];
                        if((level+1) <= totalLevel)
                        {
                            unlockRow[consolidated_ctr + 1] = levelPtr[level+1];
                        }
                        consolidated_ctr += 1;
                        consolidated_curLevelCtr += 1;
                    }
                }
            }
            curLevelCtr++;
        }
        int lastLevel = levelGroupPtr[levelGroup+1];
        newLevelPtr.push_back(levelPtr[lastLevel]);

        dangerRow[consolidated_ctr] = levelPtr[lastLevel-1];
        if((lastLevel+1) <= totalLevel)
        {
            unlockRow[consolidated_ctr + 1] = levelPtr[lastLevel+1];
        }

        consolidated_ctr += 1;
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

void mtxPower::getStatNUMA()
{
    printf("NUMA statistics\n");
    for(int i=0; i<(numSharedCache); ++i)
    {
        int NUMAnnz = 0;
        int NUMAnrow = 0;
        for(int level=levelGroupPtr[i]; level<levelGroupPtr[i+1]; ++level)
        {
            NUMAnnz += levelData->levelNnz[level];
            NUMAnrow += levelData->levelRow[level];
        }
        printf("levelGroup: %d, Level = [%d, %d], Nrow = %d, Nnz = %d\n", i, levelGroupPtr[i], levelGroupPtr[i+1], NUMAnrow, NUMAnnz);
    }
}

void mtxPower::findUnlockCtr()
{

    if(!unlockCtr)
    {
        unlockCtr = (int*)malloc(sizeof(int)*totalLevel);

        for(int l=0; l<totalLevel; ++l)
        {
            unlockCtr[l] = 0;
        }


#pragma omp parallel
        {

            int totalLevelGroup = getTotalLevelGroup();
            int threadPerLevelGroup = omp_get_num_threads()/totalLevelGroup;
            int tid = omp_get_thread_num();
            int levelGroup = tid / threadPerLevelGroup;
            int localTid = tid % threadPerLevelGroup;
            int startLevel = levelGroupPtr[levelGroup];
            int endLevel = levelGroupPtr[levelGroup+1];
            int offset = 0; //doesn't matter what is offset
            //initialize lockTableCtr
            for(int l=startLevel; l<endLevel; ++l)
            {
                SPLIT_LEVEL_PER_THREAD_P2P(l);
                if(currUnlockRow > startRow_tid)
                {
#pragma omp atomic
                    ++unlockCtr[l];
                }
                if(0)
                {
                    ERROR_PRINT("Should never be here");
                    //Suppress unused warning
                    printf("%d, %d\n", endRow_tid, dangerRowStart);
                }
            }
        }
    }
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

int* mtxPower::getUnlockRowRef()
{
    return unlockRow;
}

int* mtxPower::getDangerRowRef()
{
    return dangerRow;
}

int* mtxPower::getUnlockCtrRef()
{
    return unlockCtr;
}


