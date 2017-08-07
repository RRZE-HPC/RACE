#include "lb.h"
#include "utility.h"
#include <cmath>
#include "macros.h"

LB::LB(int nThreads_, LevelData* levelData_, dist_t dist_, d2Method d2Type_, LB_t lbTarget_):levelPtr(NULL),zonePtr(NULL),levelData(levelData_), dist(dist_), d2Type(d2Type_), maxThreads(-1),nThreads(nThreads_), lbTarget(lbTarget_)
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
    //First do NAIVE partititon
    int blockWidth = 0;
    int totalLevel = levelData->totalLevel;

    maxThreads = getPossibleThreads(totalLevel, dist, d2Type);
    blockPerThread = getBlockPerThread(dist, d2Type);
    minGap = getMinGap(dist, d2Type);

    if(nThreads > maxThreads)
    {
        WARNING_PRINT("Requested threads(%d) cannot be used, limit to %d threads",nThreads,maxThreads);
        nThreads = maxThreads;
    }

    blockWidth = int(totalLevel/(blockPerThread*nThreads));//2 blocks for 1 thread of this width
    nBlocks = blockPerThread*nThreads;
    int levelPtrSize = nBlocks+1;
    levelPtr = new int[levelPtrSize];

    levelPtr[0] = 0;
    for(int i=0; i<levelPtrSize-1; ++i) {
        levelPtr[i+1] = levelPtr[i] + blockWidth;
    }
    levelPtr[levelPtrSize-1] = totalLevel;

    //Now do load balancing
    if(nThreads < maxThreads)
    {
        int *chunkSum = new int[nBlocks];
        int *newChunkSum = new int[nBlocks];

        bool exitFlag = false;
        int *oldLevelPtr = new int[nBlocks+1];

        while(!exitFlag)
        {
            calcChunkSum(chunkSum);
            Stat meanVar(chunkSum, nBlocks, blockPerThread);
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
                    Stat newMeanVar(newChunkSum, nBlocks, blockPerThread);
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

Stat::Stat(int* arr_, int len_, int numPartitions_):arr(arr_),len(len_),numPartitions(numPartitions_),mean(NULL),var(NULL),acquireWeight(NULL),giveWeight(NULL),weight(NULL)
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
        for(int i=partition; i<len; i+=numPartitions) {
            mean[partition] = mean[partition] + arr[i];
        }
        mean[partition] = (numPartitions*mean[partition])/(len);

        for(int i=partition; i<len; i+=numPartitions) {
            double diff = (arr[i]-mean[partition]);
            var[partition] = var[partition] + diff*diff;
            acquireWeight[i] = diff;
            giveWeight[i] = -diff;
            weight[i] = std::fabs(diff);
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
