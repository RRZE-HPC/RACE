#include "matrixPower.h"
#include "utility.h"
#include <cmath>
#include "macros.h"
#include "lb.h"
#include "omp.h"

//#define LB_REMINDER

mtxPower::mtxPower(Graph* graph_, int highestPower_, int numSharedCache_, double cacheSize_, double safetyFactor_, int cache_violation_cutoff_, int startRow_, int endRow_):graph(graph_), cacheLevelGroup(NULL),levelData(NULL), highestPower(highestPower_), numSharedCache(numSharedCache_), cacheSize(cacheSize_), safetyFactor(safetyFactor_), cache_violation_cutoff(cache_violation_cutoff_), startRow(startRow_), endRow(endRow_)
{
    traverser = new Traverse(graph, RACE::ONE, startRow, endRow);
}

mtxPower::~mtxPower()
{
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
    levelPtr =  std::vector<int>(totalLevel+1);

    levelPtr[0] = startRow;

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

int mtxPower::workingBoundaryLength()
{
    return static_cast<int>(ceil(highestPower/2.0))-1;
}

void mtxPower::identifyHopelessRegions(std::vector<int> cacheViolatedLevel)
{
    if(cacheViolatedLevel.size() > 1)
    {
        cacheViolatedLevel.push_back(cacheViolatedLevel.back()+2); //a dummy so last element will be pushed, the dummy is made 2 apart from last so it doesn't compress
    }
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
                hopelessRegions.push_back(cacheViolatedLevel[i]);
            }
            else if((currFlag > 0) && (prevFlag > 0))
            {
                hopelessRegions.push_back(cacheViolatedLevel[i-1]);
                hopelessRegions.push_back(cacheViolatedLevel[i]);
            }
            else if(currFlag == 0)
            {
                //this is end of a hopeless region
                hopelessRegions.push_back(cacheViolatedLevel[i-1]);
            }

            //printf("@@@@ check hopelessRegions[%d] = %d\n", (int)hopelessRegions.size()-1, hopelessRegions.back());
        }
        prevFlag = currFlag;
    }

    for(int i=0; i<(int)hopelessRegions.size(); ++i)
    {
        printf("hopelessRegions[%d] = %d\n", i, hopelessRegions[i]);
    }
    printf("merging regions that do not have min. distance\n");
    //combine hopeless region that do not have a distance of atleast
    //(highestPower-1) between them; Note here distance of atleast totalBoundary
    //is required (and not just workingBoundary)
    int n_hopeless =  (int)(hopelessRegions.size()/2.0);
    std::vector<int> new_hopelessRegions;

    int curStart, curEnd;
    int nxtStart, nxtEnd;
    getHopelessStartEnd(0, &curStart, &curEnd);
    getHopelessStartEnd(1, &nxtStart, &nxtEnd);

    for(int i=0; i<(n_hopeless); ++i)
    {
        if(curStart != -1)
        {
            while(((nxtStart-curEnd) < (highestPower-1)) && (nxtStart!=-1))
            {
                //distance is less, merging
                curEnd = nxtEnd;
                getHopelessStartEnd(i+2, &nxtStart, &nxtEnd);
                i++;
            }
            new_hopelessRegions.push_back(curStart);
            new_hopelessRegions.push_back(curEnd);
            getHopelessStartEnd(i+1, &curStart, &curEnd);
            getHopelessStartEnd(i+2, &nxtStart, &nxtEnd);
        }
    }

    hopelessRegions = new_hopelessRegions;
    n_hopeless =  (int)(hopelessRegions.size()/2.0);

    new_hopelessRegions.resize(0);

    hopelessNodePtr.resize(numSharedCache+1,0);

    for(int i=0; i<(int)hopelessRegions.size(); ++i)
    {
        printf("hopelessRegions[%d] = %d\n", i, hopelessRegions[i]);
    }

    std::vector<int> toDrop;
    printf("split hopeless with node boundaries, and ensure min. distance from node boundaries\n");
    //make sure hopelessRegions doesn't start from nodeBoundaries,
    //and split if spawning across boundaries
    for(int node=0; node<numSharedCache; ++node)
    {
        int levelStart = nodePtr[node];
        int levelEnd = nodePtr[node+1];

        //first split the ones accross nodeBoundaries
        for(int i=0; i<n_hopeless; ++i)
        {
            getHopelessStartEnd(i, &curStart, &curEnd);
            if(curStart >= levelStart)
            {
                if(curStart < levelEnd)
                {
                //tail crossing boundary,
                    if(curEnd > levelEnd)
                    {
                        //spwaning multiple nodeBoundaries
                        new_hopelessRegions.push_back(curStart);
                        new_hopelessRegions.push_back(levelEnd-1);
                    }
                    else
                    {
                        new_hopelessRegions.push_back(curStart);
                        new_hopelessRegions.push_back(curEnd);
                    }
                }
            }
            else
            {
                if(curEnd > levelEnd)
                {
                    WARNING_PRINT("At least one complete node is hopeless\n");
                    new_hopelessRegions.push_back(levelStart);
                    new_hopelessRegions.push_back(levelEnd-1);
                }
                else
                {
                    if(curEnd >= levelStart)
                    {
                        //head crossing boundary
                        new_hopelessRegions.push_back(levelStart);
                        new_hopelessRegions.push_back(curEnd);
                    }
                }
            }
        }
        int startHopelessNodeId = hopelessNodePtr[node];
        //check if boundary within (highestPower-1), if not adjust
        getHopelessStartEnd(startHopelessNodeId, &curStart, &curEnd, new_hopelessRegions);
        //check start boundary
        if((curStart!=-1) && ((curStart-levelStart) < (highestPower-1)))
        {
            int newCurStart = levelStart+(highestPower-1);
            //Remember curStart can be equal to curEnd
            if(curEnd < newCurStart)
            {
                //drop this hopelessRegion
                //new_hopelessRegions.erase(new_hopelessRegions.begin()+2*startHopelessNodeId, new_hopelessRegions.begin()+2*startHopelessNodeId+1);
                toDrop.push_back(startHopelessNodeId);
            }
            else
            {
                new_hopelessRegions[2*startHopelessNodeId] = newCurStart;
            }
        }
        int endHopelessNodeId = static_cast<int>(new_hopelessRegions.size()/2.0)-1;
        getHopelessStartEnd(endHopelessNodeId, &curStart, &curEnd, new_hopelessRegions);
        //check end boundary; curEnd+1 since hopelessRegions include last level,
        //while levelPtr is one past last level
        if((curEnd!=-1) && ((levelEnd-(curEnd+1)) < (highestPower-1)))
        {
            int newCurEnd = levelEnd-(highestPower-1)-1;

            if(curStart > newCurEnd)
            {
                //drop this hopelessRegion
                //new_hopelessRegions.erase(new_hopelessRegions.begin()+2*endHopelessNodeId, new_hopelessRegions.begin()+2*endHopelessNodeId+1);
                toDrop.push_back(endHopelessNodeId);
            }
            else
            {
                new_hopelessRegions[2*endHopelessNodeId+1] = newCurEnd;
            }
        }

        //now find out hopelessNodePtr
        hopelessNodePtr[node+1] = static_cast<int>(new_hopelessRegions.size()/2.0);
    }
    hopelessRegions = new_hopelessRegions;
    new_hopelessRegions.resize(0);
    n_hopeless =  (int)(hopelessRegions.size()/2.0);
    std::vector<int> new_hopelessNodePtr(numSharedCache+1,0);

    for(int i=0; i<(int)hopelessRegions.size(); ++i)
    {
        printf("hopelessRegions[%d] = %d\n", i, hopelessRegions[i]);
    }

    printf("\n");
    for(int i=0; i<(int)hopelessNodePtr.size(); ++i)
    {
        printf("HopelessNodePtr[%d] = %d\n", i, hopelessNodePtr[i]);
    }


    printf("Dropping unnecessary regions\n");
    for(int i=0; i<(int)toDrop.size(); ++i)
    {
        printf("toDrop[%d] = %d\n", i, toDrop[i]);
    }

    for(int n=0; n<numSharedCache; ++n)
    {
        int startHopelessNode = hopelessNodePtr[n];
        int endHopelessNode = hopelessNodePtr[n+1];
        for(int i=startHopelessNode; i<endHopelessNode; ++i)
        {
            bool dropFlag = false;
            //now drop unwanted regions; by checking if the index is in toDrop list
            for(int j=0; j<(int)toDrop.size(); ++j)
            {
                if(i==toDrop[j])
                {
                    dropFlag = true;
                }
            }
            if(!dropFlag)
            {
                new_hopelessRegions.push_back(hopelessRegions[2*i]);
                new_hopelessRegions.push_back(hopelessRegions[2*i+1]);
            }
        }
        new_hopelessNodePtr[n+1] = static_cast<int>(new_hopelessRegions.size()/2.0);
    }
    hopelessRegions = new_hopelessRegions;
    n_hopeless =  (int)(hopelessRegions.size()/2.0);
    hopelessNodePtr = new_hopelessNodePtr;

    for(int i=0; i<(int)hopelessRegions.size(); ++i)
    {
        printf("hopelessRegions[%d] = %d\n", i, hopelessRegions[i]);
    }

    printf("\n");
    for(int i=0; i<(int)hopelessNodePtr.size(); ++i)
    {
        printf("HopelessNodePtr[%d] = %d\n", i, hopelessNodePtr[i]);
    }


    printf("save boundaries\n");
    //save the boundaries of hopeless regions
    //there are 2 kinds of boundaries
    //1) total boundary and 2) working boundary
    //1 -   total boundary includes all (highestPower-1) boundaries of the
    //      hopelessRegion where computations till highestPower is not yet reached
    //2 -   working boundary includes only the boundary region which has to be
    //      calculated with hopelessRegions. Its length is ceil(highestPower/2)-1.
    //(total boundary - working boundary)regions are done after
    //hopelessRegion is finished.
    //Only the working boundary is stored here, since the others are not
    //required to be passed to Traverse class, and while doing compuatation can
    //be obtained from levelPtr. Also others will change as consolidatePartition is called.
    hopelessRegionPositiveBoundary.resize(n_hopeless);
    hopelessRegionNegativeBoundary.resize(n_hopeless);

    int boundaryLength = workingBoundaryLength();
    for(int i=0; i<n_hopeless; ++i)
    {
        std::vector<int> curPositiveBoundary;
        std::vector<int> curNegativeBoundary;
        int curHopelessStartLevel, curHopelessEndLevel;
        getHopelessStartEnd(i, &curHopelessStartLevel, &curHopelessEndLevel);
        printf("Hopeless = [%d, %d]\n", curHopelessStartLevel, curHopelessEndLevel);
        for(int p=0; p<=boundaryLength; ++p)
        {
            int curNegativeLevel = curHopelessStartLevel-p;
            int curPositiveLevel = curHopelessEndLevel+p;

            if(curNegativeLevel < 0)
            {
                ERROR_PRINT("This shouldn't happen, hopelessRegions shouldn't be this close to boundary");
                curNegativeBoundary.push_back(-1);
            }
            else
            {
                curNegativeBoundary.push_back(levelPtr[curNegativeLevel]);
            }
            if(curPositiveLevel > (totalLevel-1))
            {
                ERROR_PRINT("This shouldn't happen, hopelessRegions shouldn't be this close to boundary");
                curPositiveBoundary.push_back(-1);
            }
            else
            {
                curPositiveBoundary.push_back(levelPtr[curPositiveLevel+1]);
            }

        }

        hopelessRegionNegativeBoundary[i] = curNegativeBoundary;
        hopelessRegionPositiveBoundary[i] = curPositiveBoundary;
    }

    //printing
    printf("Hopeless: Negative Boundary\n");
    for(int i=0; i<n_hopeless; ++i)
    {
        for(int p=0; p<=boundaryLength; ++p)
        {
            printf("%d, ", hopelessRegionNegativeBoundary[i][p]);
        }
        printf("\n");
    }
    printf("Hopeless: Positive Boundary\n");
    for(int i=0; i<n_hopeless; ++i)
    {
        for(int p=0; p<=boundaryLength; ++p)
        {
            printf("%d, ", hopelessRegionPositiveBoundary[i][p]);
        }
        printf("\n");
    }
}

//TODO: Write a wrapper function that calls findPartition
//based on levelTree info

//return levelTree
//TODO : add 2 args to findPartition: start and end
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
            if(factor > cache_violation_cutoff)
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
    identifyHopelessRegions(cacheViolatedLevel);
    consolidatePartition();
    //TODO: add the partitions to level tree
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
    lb.getZonePtr(&cacheLevelGroup, &len, startRow);

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

   // nodePtr = std::vector<int>(len);
    nodePtr = findMacroLevelPtr(cacheLevelGroup);
    if(nodePtr.size() != numSharedCache+1)
    {
        ERROR_PRINT("nodePtr dimensions do not match, nodePtr.size() = %d, numSharedCache = %d", (int)nodePtr.size(), numSharedCache);
    }

    for(int i=0; i<len; ++i)
    {
        printf("nodePtr = %d\n", nodePtr[i]);
    }

    getStatNUMA();
}

std::vector<int> mtxPower::findMacroLevelPtr(int* zones)
{
    std::vector<int> macroLevelPtr;
    int ctr = 0;
    printf("levelPtr size = %d\n", (int)levelPtr.size());
    for(int i=0; i<(totalLevel+1); ++i)
    {
        if((ctr < (numSharedCache+1)) && (levelPtr[i] == zones[ctr])) //ctr can be greater because of level with 0 rows, so guard it
        {
            macroLevelPtr.push_back(i);
            ++ctr;
        }
    }
    return macroLevelPtr;
}

void mtxPower::getHopelessStartEnd(int count, int *start, int *end)
{
    getHopelessStartEnd(count, start, end, hopelessRegions);
}

void mtxPower::getHopelessStartEnd(int count, int *start, int *end, std::vector<int> _hopelessRegions_)
{
    (*start) = -1;
    (*end) = -1;
    if((count >= 0) && (2*count <  (int)_hopelessRegions_.size()-1))
    {
        (*start) = _hopelessRegions_[2*count];
        (*end) = _hopelessRegions_[2*count+1];
    }
}

void mtxPower::consolidatePartition()
{
    unlockRow = std::vector<int>(totalLevel);
    dangerRow = std::vector<int>(totalLevel);

    std::vector<int> newLevelPtr;
    newLevelPtr.push_back(startRow); //done in loop
    int consolidated_ctr = 0; //startRow;
    unlockRow[0] = levelPtr[1];

#ifdef LB_REMINDER
    //find max start and end boundary levels
    std::vector<int> maxSumBoundaries(highestPower-1,0); //counter for each boundary level

    for(int pow=0; pow<highestPower-1; ++pow)
    {
        for(int node=0; node<numSharedCache; ++node)
        {
            int startBoundary = nodePtr[node]+pow;
            int endBoundary = nodePtr[node+1]-1-pow;
            int startNNZ = levelData->levelNnz[startBoundary];
            int endNNZ =  levelData->levelNnz[endBoundary];

            maxSumBoundaries[pow] = std::max(maxSumBoundaries[pow],startNNZ+endNNZ);
        }
        printf("Boundary @ power %d = %d\n", pow+1, maxSumBoundaries[pow]);
    }
#endif

    std::vector<int> consolidated_hopelessRegions;
    int hopelessRegionCount = 0;
    int currHopelessStart, currHopelessEnd;
    getHopelessStartEnd(hopelessRegionCount, &currHopelessStart, &currHopelessEnd);

    printf("check hopeless = [%d, %d]\n", currHopelessStart, currHopelessEnd);
    for(int node=0; node<numSharedCache; ++node)
    {
        double sumElem = 0;
        double sumNNZ = 0;
        int curLevelCtr = 0;
        int consolidated_curLevelCtr = 0;
        //int totalLevelInGroup = nodePtr[node+1]-nodePtr[node];


        //newLevelPtr.push_back(levelPtr[nodePtr[node]]);
        for(int level=nodePtr[node]; level<nodePtr[node+1]; ++level)
        {
            double currElem = getElemUpperLimit(level);
            double bytePerNNZ = getBytePerNNZ();
            double cacheElem = cacheSize/bytePerNNZ;
            sumElem += currElem;

            double curNNZ = levelData->levelNnz[level];
            sumNNZ += curNNZ;

            bool hopeless = false;
            bool forceUpdate = false;

            int boundaryLength = workingBoundaryLength();

            if(level == currHopelessStart)
            {
                consolidated_hopelessRegions.push_back(consolidated_ctr+1);//+1 since not yet incrmented
                //Only start is needed since we know it would always be one
                //level long after consolidation
                //consolidated_hopelessRegions.push_back(consolidated_ctr+2);//+1 since it will be consecuive, note that this is pointing towards end and does not include end, contrary to hopelessRegions
            }



            //check if level in working boundary of hopelessRegions; 
            //start boundary
            if( (level >= (currHopelessStart-boundaryLength)) && (level <= currHopelessStart) )
            {
                forceUpdate = true;
            }
            //end boundary
            else if((level > currHopelessEnd) && (level <= (currHopelessEnd+boundaryLength)) )
            {
                forceUpdate = true;
            }
            //deal with hopeless region
            else if((level > currHopelessStart) && (level <= currHopelessEnd))
            {
                //this level belongs to hopeless region
                hopeless = true;
            }
            /*else if(level == currHopelessEnd)
            {
                //don't set hopeless because now we have to end
                //and consolidate
                hopelessRegionCount++;
                getHopelessStartEnd(hopelessRegionCount, &currHopelessStart, &currHopelessEnd);
            }*/
            if(level == (currHopelessEnd+boundaryLength+1))
            {
                //don't set hopeless because now we have to end
                //and consolidate
                hopelessRegionCount++;
                getHopelessStartEnd(hopelessRegionCount, &currHopelessStart, &currHopelessEnd);
                printf("check hopeless = [%d, %d]\n", currHopelessStart, currHopelessEnd);
                forceUpdate = true;
            }
            //printf("#### check level=%d, forceUpdate = %d, hopeless = %d\n", level, forceUpdate, hopeless);
            //printf("##### hopeless = %d, hopeless = [%d,%d]\n", hopeless, currHopelessStart, currHopelessEnd);
#ifdef LB_REMINDER
            int startPow = 0;
#endif
            if(curLevelCtr > 0 && !hopeless)
            {
                //second condition to avoid bulk levels at boundary,
                //as this is processed without caching

#ifdef LB_REMINDER
                if( (sumElem >= cacheElem) || ( (numSharedCache > 1) && ( (consolidated_curLevelCtr < (highestPower)) || (curLevelCtr > (totalLevelInGroup-highestPower)) ) ) )
#else
                if((sumElem >= cacheElem) || forceUpdate)
#endif
                {
                    bool updateFlag = true;

#ifdef LB_REMINDER
                    int endBoundary = nodePtr[node+1]-(startPow)-1;
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
        int lastLevel = nodePtr[node+1];
        newLevelPtr.push_back(levelPtr[lastLevel]);

        dangerRow[consolidated_ctr] = levelPtr[lastLevel-1];
        if((lastLevel+1) <= totalLevel)
        {
            unlockRow[consolidated_ctr + 1] = levelPtr[lastLevel+1];
        }

        consolidated_ctr += 1;
    }

    hopelessRegions = consolidated_hopelessRegions;

    for(int i=0; i<consolidated_hopelessRegions.size(); ++i)
    {
        printf("consolidated_hopelessRegions[%d] = %d\n", i, consolidated_hopelessRegions[i]);
    }

    //rewrite levelPtr with newLevelPtr
    totalLevel = (int)(newLevelPtr.size()-1);
    levelPtr.resize(totalLevel+1);
    printf("consolidated levels = %d\n", totalLevel);
    for(int i=0; i<totalLevel+1; ++i)
    {
        levelPtr[i] = newLevelPtr[i];
        printf("levelPtr[%d] = %d\n", i, levelPtr[i]);
    }

    std::vector<int> newNodePtr = findMacroLevelPtr(cacheLevelGroup);
    if(nodePtr.size() != numSharedCache+1)
    {
        ERROR_PRINT("nodePtr dimensions do not match");
    }


    for(int i=0; i<(numSharedCache+1); ++i)
    {
        nodePtr[i] = newNodePtr[i];
        printf("nodePtr_consolidated[%d] = %d, row = %d\n", i, nodePtr[i], levelPtr[nodePtr[i]]);
    }

}

void mtxPower::getStatNUMA()
{
    printf("NUMA statistics\n");
    for(int i=0; i<(numSharedCache); ++i)
    {
        int NUMAnnz = 0;
        int NUMAnrow = 0;
        for(int level=nodePtr[i]; level<nodePtr[i+1]; ++level)
        {
            NUMAnnz += levelData->levelNnz[level];
            NUMAnrow += levelData->levelRow[level];
        }
        printf("node: %d, Level = [%d, %d], Nrow = %d, Nnz = %d\n", i, nodePtr[i], nodePtr[i+1], NUMAnrow, NUMAnnz);
    }
}

void mtxPower::findUnlockCtr()
{

    if(unlockCtr.empty())
    {
        unlockCtr = std::vector<int>(totalLevel, 0);

#pragma omp parallel
        {

            int totalNodes = getTotalNodes();
            int threadPerNode = omp_get_num_threads()/totalNodes;
            int tid = omp_get_thread_num();
            int node = tid / threadPerNode;
            int localTid = tid % threadPerNode;
            int startLevel = nodePtr[node];
            int endLevel = nodePtr[node+1];
            int offset = 0; //doesn't matter what is offset
            //initialize lockTableCtr
            for(int l=startLevel; l<endLevel; ++l)
            {
                SPLIT_LEVEL_PER_THREAD_P2P_NOREF(l);
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

std::vector<int> mtxPower::getLevelPtr()
{
    return levelPtr;
}

int mtxPower::getTotalNodes()
{
    return numSharedCache;
}

std::vector<int> mtxPower::getNodePtr()
{
    return nodePtr;
}

std::vector<int> mtxPower::getUnlockRow()
{
    return unlockRow;
}

std::vector<int> mtxPower::getDangerRow()
{
    return dangerRow;
}

std::vector<int> mtxPower::getUnlockCtr()
{
    return unlockCtr;
}

std::vector<int> mtxPower::getHopelessRegions()
{
    return hopelessRegions;
}

std::vector<int> mtxPower::getHopelessNodePtr()
{
    return hopelessNodePtr;
}
