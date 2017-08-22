#include "lb.h"
#include "utility.h"
#include <cmath>
#include "macros.h"

LB::LB(int nThreads_, double efficiency_, LevelData* levelData_, dist_t dist_, d2Method d2Type_, LB_t lbTarget_):levelPtr(NULL),zonePtr(NULL),levelData(levelData_), dist(dist_), d2Type(d2Type_), maxThreads(-1), nThreads(nThreads_), efficiency(efficiency_), lbTarget(lbTarget_)
{
    printf("dist = %d\n",dist);
    if( (lbTarget == NNZ) && (levelData->levelNnz == NULL) )
    {
        WARNING_PRINT("levelDataNnz not recieved, Load balancing target will fallback to ROW");
        lbTarget = ROW;
    }
}

LB::~LB()
{
    if(levelPtr)
    {
        delete[] levelPtr;
    }

    if(zonePtr)
    {
        delete[] zonePtr;
    }

    if(effRatio)
    {
        delete[] effRatio;
    }

}

void LB::calcChunkSum(int *arr, bool forceRow)
{
    int *targetLevelData = NULL;
    if( (lbTarget == ROW) || forceRow) {
        targetLevelData = levelData->levelRow;
    } else {
        targetLevelData = levelData->levelNnz;
    }

    for(int i=0; i<nBlocks; ++i)
    {
        arr[i] = 0;
        for(int j=levelPtr[i]; j<levelPtr[i+1]; ++j)
        {
            arr[i] += targetLevelData[j];
        }
    }
}

void LB::calcEffectiveRatio()
{
    effRatio = new double[levelData->totalLevel];
    int *targetLevelData = NULL;
    int targetNElements = 0;
    if(lbTarget == ROW) {
        targetLevelData = levelData->levelRow;
        targetNElements = levelData->nrow;
    } else {
        targetLevelData = levelData->levelNnz;
        targetNElements = levelData->nnz;
    }

    for(int i=0; i<levelData->totalLevel; ++i)
    {
        effRatio[i] = (targetLevelData[i]*nThreads)/((double)targetNElements);
    }
}


int LB::findNeighbour(const Stat &stats, neighbour_t type)
{
    int *perm = new int[nBlocks];
    int neighbourIdx = -1;

    for(int i=0; i<nBlocks; ++i)
    {
        perm[i] = i;
    }

    if(type == acquire) {
        sortPerm(stats.acquireWeight, perm, 0, nBlocks, true);
    } else {
        sortPerm(stats.giveWeight, perm, 0, nBlocks, true);
    }

    int rankIdx = 0;

    if(type == acquire)
    {
        while(rankIdx<nBlocks)
        {
            int currIdx = perm[rankIdx];
            if((levelPtr[currIdx+1] - levelPtr[currIdx]) > minGap) {
                neighbourIdx = currIdx;
                break;
            }
            ++rankIdx;
        }
    } else {
        neighbourIdx = perm[0];
    }

    delete[] perm;
    return neighbourIdx;
}

void LB::moveOneStep(int toIdx, int fromIdx)
{
    if(toIdx<fromIdx)
    {
        for(int i=toIdx+1; i<=fromIdx; ++i)
        {
            levelPtr[i] += 1;
        }
    } else {
        for(int i=toIdx; i>=fromIdx+1; --i)
        {
            levelPtr[i] -= 1;
        }
    }
}


void LB::splitZones()
{
    int totalLevel = levelData->totalLevel;

/*    maxThreads = getPossibleThreads(totalLevel, dist, d2Type);
    blockPerThread = getBlockPerThread(dist, d2Type);
    minGap = getMinGap(dist, d2Type);

    if(nThreads > maxThreads)
    {
        WARNING_PRINT("Requested threads(%d) cannot be used, limit to %d threads",nThreads,maxThreads);
        nThreads = maxThreads;
    }
*/
    calcEffectiveRatio();

    minGap = getMinGap(dist, d2Type);
    blockPerThread = getBlockPerThread(dist, d2Type);

    int minGapPerThread = minGap*blockPerThread;
    double accumulator = 0;
    int minGapCtr = 0;

    double currRem = 0;
    double localRem = 0;
    double newLocalRem = 0;
    double minLocalRem = 0;
    int bestIdx = 0;
    std::vector<int> initialLvlPtr;
    initialLvlPtr.push_back(0);

    //calculate scale values
    for(int i=0; i<totalLevel; ++i)
    {
        accumulator += effRatio[i];
        ++minGapCtr;

        if(minGapCtr >= minGapPerThread)
        {
            if(roundDouble(accumulator) != 0)
            {
                if( ((1-std::fabs(accumulator-roundDouble(accumulator))) *100) > efficiency )
                {
                    //the criteria is accepted; but we can check to get more
                    //efficiency and if possible  to reduce reminder
                    localRem =  (accumulator-roundDouble(accumulator));
                    newLocalRem = localRem;
                    minLocalRem = localRem;
                    int roundedVal = roundDouble(accumulator);

                    double tempAccumulator = accumulator;
                    double minTempAccumulator = accumulator;
                    int j = i+1;

                    bestIdx = i;
                    //check for other nearest point close to rounded value
                    while( (j<totalLevel) && (std::fabs(newLocalRem)<=std::fabs(localRem)) )
                    {
                        tempAccumulator += effRatio[j];
                        newLocalRem = (tempAccumulator-roundedVal);
                        if(std::fabs(minLocalRem) > std::fabs(newLocalRem))
                        {
                            minLocalRem = newLocalRem;
                            minTempAccumulator = tempAccumulator;
                            bestIdx = j;
                        }
                        //if sufficiently close
                        else if( (std::fabs(minLocalRem) - std::fabs(newLocalRem)) < 1e-6)
                        {
                            double lhsCurrRem = minLocalRem + currRem;
                            double rhsCurrRem = newLocalRem + currRem;

                            if(lhsCurrRem > rhsCurrRem)
                            {
                                minLocalRem = newLocalRem;
                                minTempAccumulator = tempAccumulator;
                                bestIdx = j;
                            }
                        }
                        ++j;
                    }

                    //retrieve bestIdx
                    minGapCtr += (bestIdx-i);
                    i = bestIdx;
                    accumulator = minTempAccumulator;

                    printf("i=%d accumulator = %f\n", i, accumulator);

                    currRem = currRem + (accumulator-roundDouble(accumulator));
                    scale.push_back(roundDouble(accumulator));
                    initialLvlPtr.push_back(i+1);
                    //to compensate for reminder
                    accumulator = 0;
                    minGapCtr = 0;
                }
            }
        }

        if(i>=(totalLevel-5))
        {
            //push everything to the last scale value
            //else last thread won't get min. levels
            double lastThreadAccumulator = accumulator;

            for(int j=i+1; j<totalLevel; ++j)
            {
                lastThreadAccumulator += effRatio[j];
            }
/*          int currThreads = 0;
            for(int k=0; k<scale.size(); ++k)
            {
                currThreads += scale[k];
            }*/
            if(roundDouble(lastThreadAccumulator)!=0)
            {
                scale.push_back(roundDouble(lastThreadAccumulator));
                initialLvlPtr.push_back(totalLevel);
            }
            else
            {
                initialLvlPtr.back() = totalLevel;
            }
            break;
        }
    }

    //do a check to ensure correct number of threads are produced
    //for initial estimate
    int threadCtr = 0;
    for(int i=0; i!=scale.size(); ++i)
    {
        if(scale[i]==0)
        {
            scale[i] = 1;
        }
        threadCtr+=scale[i];
    }

    if(threadCtr != nThreads)
    {
        WARNING_PRINT("This message will be switched off: Diff. in Initial estimate spawned threads = %d, required = %d",threadCtr,nThreads );
        int diff = nThreads - threadCtr;
        //adjust this difference; this might 
        //have occured due to rounding
        if(std::abs(diff)==1)
        {
            if(diff==1)
            {
                scale[scale.size()-1] += 1;
            }
            else
            {
                if(scale[scale.size()-1] == 1)
                {
                    scale.pop_back();
                    initialLvlPtr.pop_back();
                    initialLvlPtr[initialLvlPtr.size()-1] = totalLevel;
                }
                else
                {
                    scale[scale.size()-1] -= 1;
                }
            }
        }
        else
        {
            WARNING_PRINT("Given number of threads would not be spawned; please reduce efficiency if more threads are needed");
        }
    }

    currLvlThreads = scale.size();

    for(int i=0; i<currLvlThreads; ++i)
    {
        printf("scale[%d] = %d\n", i, scale[i]);
    }

    for(int i=0; i<currLvlThreads+1; ++i)
    {
        printf("initLvlPtr[%d] = %d\n", i, initialLvlPtr[i]);
    }


    nBlocks = blockPerThread*currLvlThreads;
    int levelPtrSize = nBlocks+1;
    levelPtr = new int[levelPtrSize];

    levelPtr[0] = 0;
    for(int i=0; i<currLvlThreads; ++i)
    {
        for(int j=1; j<=blockPerThread; ++j)
        {
            levelPtr[i*blockPerThread+j] = initialLvlPtr[i] + roundDouble(j*(initialLvlPtr[i+1]-initialLvlPtr[i])/((double)blockPerThread));
        }
    }

    for(int i=0; i<levelPtrSize; ++i)
    {
        printf("levelPtr[%d] = %d\n",i,levelPtr[i]);
    }

    //Now do load balancing; with scale values
    //therefore we reserve for nested threads
    //here itself
    int *chunkSum = new int[nBlocks];
    int *newChunkSum = new int[nBlocks];

    bool exitFlag = false;
    int *oldLevelPtr = new int[nBlocks+1];

    while(!exitFlag)
    {
        calcChunkSum(chunkSum);
        Stat meanVar(chunkSum, nBlocks, blockPerThread, scale);
        double var = 0;
        for(int i=0; i<meanVar.numPartitions; ++i)
        {
            var += meanVar.var[i];
        }
        double newVar = var;

        int *rankPerm = new int[nBlocks];
        for(int i=0; i<nBlocks; ++i)
        {
            rankPerm[i] = i;
        }

        sortPerm(meanVar.weight, rankPerm, 0, nBlocks, true);
        int currRank = 0;

        for(int i=0; i<nBlocks+1; ++i)
        {
            oldLevelPtr[i] = levelPtr[i];
        }

        while(newVar>=var)
        {
            int rankIdx = rankPerm[currRank];
            for(int i=0; i<nBlocks+1; ++i)
            {
                levelPtr[i] = oldLevelPtr[i];
            }
            //determine the mean of the rank: TODO: for 3 block case
            double myMean = meanVar.mean[rankIdx%blockPerThread];
            bool fail = false;

            //If I am less than mean
            if(chunkSum[rankIdx] < myMean)
            {
                //try to acquire from my neighbours
                int acquireIdx = findNeighbour(meanVar, acquire);
                if(acquireIdx==-1)
                {
                    fail = true;
                }
                moveOneStep(rankIdx, acquireIdx);
            }
            // If I am greater than mean
            else if( (levelPtr[rankIdx+1] - levelPtr[rankIdx]) > minGap)
            {
                //try to give to my neighbours
                int giveIdx = findNeighbour(meanVar, give);
                if(giveIdx==-1)
                {
                    fail = true;
                }
                moveOneStep(giveIdx, rankIdx);
            }

            if(!fail)
            {
                calcChunkSum(newChunkSum);
                Stat newMeanVar(newChunkSum, nBlocks, blockPerThread, scale);
                //newMeanVar.calculate();

                newVar = 0;
                for(int i=0; i<newMeanVar.numPartitions; ++i)
                {
                    newVar += newMeanVar.var[i];
                }
            }

            if( (currRank == (nBlocks-1)) && (newVar>=var) )
            {
                exitFlag = true;
                delete[] levelPtr;
                levelPtr = oldLevelPtr;
                break;
            }
            currRank += 1;
        }
        delete[] rankPerm;
    }
    delete[] chunkSum;
    delete[] newChunkSum;


    printf("Final lvl ptr\n");
    for(int i=0; i<levelPtrSize; ++i)
    {
        printf("levelPtr[%d] = %d\n",i,levelPtr[i]);
    }

    calcZonePtr(0);

    for(int i=0; i<nBlocks; ++i)
    {
        printf("%d (%f)\n", (zonePtr[i+1]-zonePtr[i]), (zonePtr[i+1]-zonePtr[i])*100.0/((double)levelData->nrow));
    }
}

void LB::calcZonePtr(int base)
{
    int *chunkSum = new int [nBlocks];
    calcChunkSum(chunkSum, true);

    zonePtr = new int [nBlocks+1];

    zonePtr[0] = base;
    for(int i=0; i<nBlocks; ++i)
    {
        zonePtr[i+1] = zonePtr[i] + chunkSum[i];
    }

    delete[] chunkSum;
}

void LB::getZonePtr(int **zonePtr_, int *len, int base)
{
    calcZonePtr(base);
    (*zonePtr_) = zonePtr;
    zonePtr = NULL;
    (*len) = nBlocks+1;
}

/*The scale value can be used
 * to get next level number of threads.
 * It's an array since each thread in current
 * Lvl can spawn different amount of threads
 * in nxtLvl.
 */
void LB::getNxtLvlThreads(int **nxtLvlThreads, int *len)
{
    (*len) = scale.size();
    (*nxtLvlThreads) = new int[(*len)];

    for(int i=0; i<(*len); ++i)
    {
        (*nxtLvlThreads)[i] = scale[i];
    }
}

/*returns thread in current lvl*/
int LB::getCurrLvlNThreads()
{
    return scale.size();
}

Stat::Stat(int* arr_, int len_, int numPartitions_, std::vector<int> scale_):arr(arr_),len(len_),numPartitions(numPartitions_),scale(scale_),mean(NULL),var(NULL),acquireWeight(NULL),giveWeight(NULL),weight(NULL)
{
    mean = new double[numPartitions];
    var = new double[numPartitions];
    weight = new double[len];
    acquireWeight = new double[len];
    giveWeight = new double[len];

    calculate();
}

//TODO error handling
NAME_error LB::balance()
{
    splitZones();
    return NAME_SUCCESS;
}

void Stat::calculate()
{
    for(int partition=0; partition<numPartitions; ++partition) {
        mean[partition] = 0;
        var[partition] = 0;
        int ctr = 0;
        for(int i=partition; i<len; i+=numPartitions) {
            mean[partition] = mean[partition] + arr[i]/((double)scale[ctr]);
            ++ctr;
        }
        mean[partition] = (numPartitions*mean[partition])/(len);

        ctr = 0;
        for(int i=partition; i<len; i+=numPartitions) {
            double diff = ( (arr[i]/((double)scale[ctr]))-mean[partition]);
            var[partition] = var[partition] + diff*diff;
            acquireWeight[i] = diff;
            giveWeight[i] = -diff;
            weight[i] = std::fabs(diff);
            ++ctr;
        }
    }
}

int LB::getNumThreads()
{
    return nThreads;
}

int LB::getMaxThreads()
{
    return maxThreads;
}


Stat::~Stat()
{
    delete[] mean;
    delete[] var;
    delete[] weight;
    delete[] acquireWeight;
    delete[] giveWeight;
}
