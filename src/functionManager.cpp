#include "functionManager.h"
#include "test.h"
#include <iostream>
#include <typeinfo>
#include "machine.h"
#include "macros.h"
#include <emmintrin.h>

#define POWER_WITH_FLUSH_LOCK
//#define RACE_DEBUG

//coloring
FuncManager::FuncManager(void (*f_) (int,int,void*), void *args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart_):rev(false),power_fn(false), numaSplit(false), func(f_),args(args_),zoneTree(zoneTree_), pool(pool_), serialPart(serialPart_)
{
    //    func = std::bind(f_,*this,std::placeholders::_1,std::placeholders::_2,args);
    /*len = 72*72*72;
      a = new double[len];
      b = new double[len];
      c = new double[len];
      d = new double[len];
      */
    //  barrierTime = 0;
    recursiveFun = std::bind(recursiveCall, this, std::placeholders::_1);
#pragma omp parallel
    {
        activethreads = omp_get_num_threads();
    }
}

//power = true, NUMA = false
FuncManager::FuncManager(void (*f_) (int,int,int,int,void*), void *args_, int power_, mtxPowerRecursive *matPower_, bool numaSplit_):power_fn(true), numaSplit(numaSplit_), powerFunc(f_), args(args_), power(power_), matPower(matPower_)
{
#ifdef POWER_WITH_FLUSH_LOCK
    //    initPowerRun();
    lockCtr = std::vector<volatile int*>(matPower->tree.size(), NULL);
    lockTableCtr = std::vector<volatile int*>(matPower->tree.size(), NULL);
#endif
#pragma omp parallel
    {
        activethreads = omp_get_num_threads();
    }

    recursivePowerFun = std::bind(recursivePowerCall, this, std::placeholders::_1);

}

FuncManager::FuncManager(const FuncManager &obj)
{
    rev = obj.rev;
    power_fn = obj.power_fn;
    func = obj.func;
    powerFunc = obj.powerFunc;
    args = obj.args;
    power = obj.power;
    zoneTree = obj.zoneTree;
    pool = obj.pool;
    serialPart = obj.serialPart;
    lockCtr = obj.lockCtr;
    lockTableCtr = obj.lockTableCtr;
    numaSplit = obj.numaSplit;
}

FuncManager::~FuncManager()
{
    /*  delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
        */
    if(power_fn)
    {
        for(int i=0; i<(int)lockCtr.size(); ++i)
        {
            if(lockCtr[i])
            {
                delete[] lockCtr[i];
            }
#ifdef POWER_WITH_FLUSH_LOCK
            else
            {
                ERROR_PRINT("Lock counter [%d] does not exist", i);
            }
#endif
        }
    }
    if(power_fn)
    {
        for(int i=0; i<(int)lockTableCtr.size(); ++i)
        {
            if(lockTableCtr[i])
            {
                free((void*)lockTableCtr[i]);
            }
#ifdef POWER_WITH_FLUSH_LOCK
            else
            {
                ERROR_PRINT("Lock table counter [%d] does not exist", i);
            }
#endif
        }
    }
}


#ifdef RACE_KERNEL_THREAD_OMP
//OMP version nested pinning not working
void recursiveCall(FuncManager* funMan, int parent)
{
    std::vector<int> *children = &(funMan->zoneTree->at(parent).children);
    int totalSubBlocks = funMan->zoneTree->at(parent).totalSubBlocks;
    int blockPerThread = getBlockPerThread(funMan->zoneTree->dist, funMan->zoneTree->d2Type);

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        int currNthreads;
        if(!children->empty())
        {
            currNthreads = (children->at(2*subBlock+1)-children->at(2*subBlock))/blockPerThread;
        }
        else
        {
            currNthreads = 0;
        }

        if(currNthreads <= 1)
        {
            std::vector<int> range = funMan->zoneTree->at(parent).valueZ;
            funMan->func(range[0],range[1],funMan->args);
        }
        else
        {
#pragma omp parallel num_threads(currNthreads)
            {
                int startBlock = 0;
                int endBlock = blockPerThread;
                int inc = 1;

                if(funMan->rev)
                {
                    startBlock = blockPerThread-1;
                    endBlock = -1;
                    inc = -1;
                }
                int tid = omp_get_thread_num();
                //Pin in each call
                //int pinOrder = zoneTree->at(children->at(2*subBlock)+2*tid).pinOrder;
                //printf("omp_proc_bind = %d\n", omp_get_proc_bind());
                //pool->pin.pinThread(pinOrder);
                for(int block=startBlock; block!=endBlock; block+=inc)
                {
                    funMan->recursiveFun(children->at(2*subBlock)+blockPerThread*tid+block);
#pragma omp barrier
                }
                //printf("Black: tid = %d, pinOrder = %d, cpu = %d\n",tid,pinOrder,sched_getcpu());
            }
        }
    }
}
#else
//PTHREAD implementation
void recursiveCall(FuncManager* funMan, int parentIdx)
{
    /*if(funMan->zoneTree->at(parentIdx).pinOrder != ((int) sched_getcpu()))
      {
      printf("Pinning Error: parentIdx=%dto %d; really pinned to cpu = %d\n", parentIdx, funMan->zoneTree->at(parentIdx).pinOrder, sched_getcpu());
      }*/

    std::vector<int>* children = &(funMan->zoneTree->at(parentIdx).children);
    int blockPerThread = getBlockPerThread(funMan->zoneTree->dist, funMan->zoneTree->d2Type);
    int totalSubBlocks = funMan->zoneTree->at(parentIdx).totalSubBlocks;

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        int nthreads;
        if(!children->empty())
        {
            nthreads = (children->at(2*subBlock+1) - children->at(2*subBlock))/blockPerThread;
        }
        else
        {
            nthreads = 0;
        }
        if(nthreads <= 1)
        {
            //        int pinOrder = zoneTree->at(parentIdx).pinOrder;
            //        printf(" pinOrder = %d, cpu = %d\n",pinOrder,sched_getcpu());
            std::vector<int>* range = &funMan->zoneTree->at(parentIdx).valueZ;
            START_TIME(func);
            funMan->func(range->front(),range->back(), funMan->args);
            STOP_TIME(func);
            funMan->zoneTree->at(parentIdx).time += GET_TIME(func);
            //test(funMan->a,funMan->b,funMan->c,funMan->d,range->at(0), range->at(1),1);
            //  test(1);
            // sleep(1);
        }
        else
        {
            int startBlock = 0;
            int endBlock = blockPerThread;
            int inc = 1;

            if(funMan->rev)
            {
                startBlock = blockPerThread-1;
                endBlock = -1;
                inc = -1;
            }
            for(int block=startBlock; block!=endBlock; block+=inc)
            {
                for(int tid=0; tid<nthreads; ++tid)
                {
                    //pool->tree[parentIdx].addJob(tid, std::bind(test,a,b,c,d, tid*len/2.0, (tid+1)*len/2.0, std::placeholders::_1), 1);
                    funMan->pool->tree[funMan->pool->poolTreeIdx(parentIdx, subBlock)].addJob(tid, funMan->recursiveFun, children->at(2*subBlock)+blockPerThread*tid+block);
                }
                funMan->pool->tree[funMan->pool->poolTreeIdx(parentIdx, subBlock)].barrier();
            }
            /*if(parentIdx == 0)
              {
              funMan->barrierTime=funMan->pool->tree[parentIdx].barrierTime;
              }*/
        }
    }
}

#endif

#if 0
//OMP nested parallelism is not a good idea for some reason when NUMA comes to play

void FuncManager::powerRun()
{
    //int totalLevel = matPower->getTotalLevel();
    int totalNodes = matPower->getTotalNodes();
    int* levelPtr = matPower->getLevelPtrRef();
    int* nodePtr = matPower->getNodePtrRef();

#pragma omp parallel num_threads(totalNodes)
    {
        int nodeGroup = omp_get_thread_num();
        //body
        {
            //if(nodeGroup == 3)
            {
                //printf("here is %d\n", sched_getcpu());
                int startLevel = nodePtr[nodeGroup];
                int endLevel = nodePtr[nodeGroup+1];
                for(int level=startLevel; level<(endLevel); ++level)
                {
                    for(int pow=0; pow<power; ++pow)
                    {
                        int powLevel = (level-pow);

                        if( (powLevel >= (startLevel+pow)) && (powLevel < (endLevel-pow)) )
                        {
                            //can be a function ptr
                            powerFunc(levelPtr[powLevel], levelPtr[powLevel+1], pow+1, args);
#if 0
#pragma omp parallel for schedule(static)
                            for(int row=levelPtr[powLevel]; row<levelPtr[powLevel+1]; ++row)
                            {
                                double tmp = 0;
                                for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
                                {
                                    tmp += A[idx]*x[(pow)*graph->NROW+col[idx]];
                                }
                                x[(pow+1)*graph->NROW+row] = tmp;
                            }
#endif
                        }
                    }
                }
            }
        }
#pragma omp barrier
        //reminder
        for(int pow=1; pow<power; ++pow)
        {
            //reminder-head
            {
                int startLevel = nodePtr[nodeGroup];
                int endLevel = std::min(startLevel+pow, nodePtr[nodeGroup+1]);
                for(int level=startLevel; level<endLevel; ++level)
                {
                    //can be a function ptr
                    powerFunc(levelPtr[level], levelPtr[level+1], pow+1, args);
                }
            }
            //reminder-tail
            {
                int endLevel = nodePtr[nodeGroup+1];
                int startLevel = std::max(nodePtr[nodeGroup], endLevel-pow);
                for(int level=startLevel; level<endLevel; ++level)
                {
                    //can be a function ptr
                    powerFunc(levelPtr[level], levelPtr[level+1], pow+1, args);
                }

            }
#pragma omp barrier
        }
    }
}

#endif

void FuncManager::initPowerRun(int parent)
{
    //need to move this routine to matrixPower.cpp
#ifdef POWER_WITH_FLUSH_LOCK
    std::vector<int>* nodePtr = &(matPower->tree[parent].nodePtr);
    int totalNodes = nodePtr->size()-1;
    int totalLevels = matPower->tree[parent].levelPtr.size()-1;

#pragma omp single
    {
        if(lockCtr[parent])
        {
            free((void*)lockCtr[parent]);
        }
        //allocate only if needed
        //if(!lockCtr[parent])
        {
            lockCtr[parent] = new volatile int[totalLevels];
        }
    }

    //TODO parallelise this routines
//#pragma omp parallel
    {
        int threadPerNode = omp_get_num_threads()/totalNodes;
        int tid = omp_get_thread_num();
        int nodeGroup = tid / threadPerNode;
        int localTid = tid % threadPerNode;
        //printf("here is %d\n", sched_getcpu());
        int startLevel = nodePtr->at(nodeGroup);
        int endLevel = nodePtr->at(nodeGroup+1);

        //only 1 thread per group, maintain NUMA
        if(localTid == 0)
        {
            for(int level=startLevel; level<endLevel; ++level)
            {
                lockCtr[parent][level] = 0;
            }
        }
    }

#pragma omp single
    {
        //pow-level semi locks
        if(lockTableCtr[parent])
        {
            free((void*)lockTableCtr[parent]);
        }
        lockTableCtr[parent] = (volatile int*)malloc(sizeof(volatile int)*power*totalLevels);

        for(int l=0; l<totalLevels; ++l)
        {
            for(int p=0; p<power; ++p)
            {
                lockTableCtr[parent][l*power+p] = 0;
            }
        }
    }

#endif
}

void FuncManager::NUMAInitPower()
{
    std::vector<int> nodePtr = matPower->tree[0].nodePtr;
    std::vector<int> levelPtr = matPower->tree[0].levelPtr;
    int totalNodes = nodePtr.size()-1;
#pragma omp parallel
    {
        int threadPerNode = omp_get_num_threads()/totalNodes;
        int tid = omp_get_thread_num();
        int nodeGroup = tid / threadPerNode;
        int numaLocalArg = (numaSplit)?nodeGroup:0;
        //only first thread in each nodeGroup works
        if(tid == nodeGroup*threadPerNode)
        {
            powerFunc(levelPtr[nodePtr[nodeGroup]], levelPtr[nodePtr[nodeGroup+1]], 1, numaLocalArg, args);
            /*for(int row=nodePtr[nodeGroup]; row<nodePtr[nodeGroup+1]; ++row)
              {
              for(int idx=rowPtr_[row]; idx<rowPtr_[row+1]; ++idx)
              {
              val[idx] = 0;
              col[idx] = 0;
              }
              }*/
        }
    }
}

//check with volatile too
#ifdef POWER_WITH_FLUSH_LOCK

#define UNLOCK(_level_, _pow_)\
{\
    if((_level_ >= 0) && (_pow_ <= power-1))\
    {\
        _Pragma("omp atomic")\
        lockTableCtr[parent][_level_*power+_pow_]++;\
    }\
}

#define WAIT(_level_, _pow_)\
{\
    if((_pow_ != 0) && (_level_!=(endLevel-1)))\
    {\
        while(lockTableCtr[parent][_level_*power+_pow_] != unlockCtr->at(_level_+1))\
        {\
            _mm_pause();\
            /*printf("Waiting on %d, need to reach%d, pow = %d, level = %d, startRow_tid = %d, unlock = %d, danger = %d, endRow_tid = %d\n", lockTableCtr[_level_*power+_pow_], unlockCtr[_level_], _pow_, _level_, startRow_tid, currUnlockRow, dangerRowStart, endRow_tid);*/\
        }\
    }\
}\


//startSlope and endSlope of computation
void FuncManager::powerCallGeneral(int startLevel, int endLevel, int boundaryStart, int boundaryEnd, int startSlope, int endSlope, int threadPerNode, const std::vector<int> *levelPtr, const std::vector<int> *unlockRow, const std::vector<int> *unlockCtr, const std::vector<int> *dangerRow, int numaLocalArg, int offset, int parent)
{
    int tid = omp_get_thread_num();
    int localTid = tid % threadPerNode;
    //printf("here is %d\n", sched_getcpu());
    /* int maxLevelCount = 0;
       for(int i=0; i<totalNodes; ++i)
       {
       maxLevelCount = std::max(maxLevelCount, (nodePtr->at(i+1)-nodePtr->at(i)));
       }
       int maxEndLevel = startLevel + maxLevelCount; //needed so that everyone calls barrier
       */

    for(int level=startLevel; level<(endLevel+power-1); ++level)
    {
        //printf("Level = %d\n", level);
        //tid responsible for (level-1) lock, locks the level
        for(int pow=0; pow<power; ++pow)
        {
            int powLevel = (level-pow);

            int boundaryLeft = boundaryStart + startSlope*pow;
            int boundaryRight = boundaryEnd + endSlope*pow;
            //wave
            if( (powLevel >= boundaryLeft) && (powLevel < boundaryRight) )
            {
                //wait till powLevel is free; use ctr[powLevel] == threads*(pow)
                while(lockCtr[parent][powLevel] < threadPerNode*pow)
                {
                    _mm_pause();
#ifdef RACE_DEBUG
                    printf("lock ctr = %d, pow = %d, level = %d\n", lockCtr[parent][powLevel], pow, powLevel);
#endif
                }

                SPLIT_LEVEL_PER_THREAD_P2P(powLevel);
#ifdef RACE_DEBUG
                printf("power = %d, pow = %d, level = %d, tid = %d, start row = %d, unlock row = %d, danger row = %d, end_row = %d\n", power, pow, powLevel, omp_get_thread_num(), startRow_tid, currUnlockRow, dangerRowStart, endRow_tid);
#endif

                //int z = (pow) + (powLevel);
                //int cons_term = std::max(0,z-power); //the term for having consequtive numbers
                //int powLevel_to_unlock = ((z*(z+1))>>1) + (pow+1) - (cons_term*(cons_term+1)>>1);

                if(dangerRowStart <= startRow_tid)
                {
                    //check lock
                    WAIT(powLevel, pow);
                    powerFunc(startRow_tid, endRow_tid, pow+1, numaLocalArg, args);
#ifdef RACE_DEBUG
                    printf("5. call tid = %d, [%d, %d], pow = %d, level = %d\n", omp_get_thread_num(), startRow_tid, endRow_tid, pow, powLevel);
#endif
                    if(startRow_tid < currUnlockRow)
                    {
                        //unlock
                        UNLOCK((powLevel-1), (pow+1));
#ifdef RACE_DEBUG
                        printf("tid = %d, Unlocking Level = %d, pow = %d\n", omp_get_thread_num(), powLevel-1, pow+1);
#endif
                    }
                }
                else if(dangerRowStart >= currUnlockRow)
                {
                    int till_row = startRow_tid;
                    if(startRow_tid < currUnlockRow)
                    {
                        till_row = std::min(currUnlockRow, endRow_tid);
                        powerFunc(startRow_tid, till_row, pow+1, numaLocalArg, args);
#ifdef RACE_DEBUG
                        printf("1. call tid = %d, [%d, %d], pow = %d, level = %d\n", omp_get_thread_num(), startRow_tid, currUnlockRow, pow, powLevel);
#endif
                        //unlock

                        UNLOCK((powLevel-1), (pow+1));
#ifdef RACE_DEBUG
                        printf("tid = %d, Unlocking Level = %d, pow = %d\n", omp_get_thread_num(), powLevel-1, pow+1);
#endif
                    }
                    /*else
                      {
                      till_row = startRow_tid;
                    //powerFunc(startRow_tid, endRow_tid, pow+1, args);
                    //till_row = endRow_tid;
                    }*/
                    if(dangerRowStart > endRow_tid)
                    {
                        powerFunc(till_row, endRow_tid, pow+1, numaLocalArg, args);
#ifdef RACE_DEBUG
                        printf("2. call tid = %d, [%d, %d], pow = %d, level = %d\n", omp_get_thread_num(), till_row, endRow_tid, pow, powLevel);
#endif
                    }
                    else
                    {
                        powerFunc(till_row, dangerRowStart, pow+1, numaLocalArg, args);
#ifdef RACE_DEBUG
                        printf("3. call tid = %d, [%d, %d], pow = %d, level = %d\n", omp_get_thread_num(), till_row, dangerRowStart, pow, powLevel);
#endif

                        //check lock
                        WAIT(powLevel, pow);

                        powerFunc(dangerRowStart, endRow_tid, pow+1, numaLocalArg, args);
#ifdef RACE_DEBUG
                        printf("4. call tid = %d, [%d, %d], pow = %d, level = %d\n", omp_get_thread_num(), dangerRowStart, endRow_tid, pow, powLevel);
#endif

                    }
                }
                else
                {
                    ERROR_PRINT("I thought this wouldn't happen, please report, tid = %d", omp_get_thread_num());
                    exit(-1);
                }
#pragma omp atomic //atomic is update by default
                lockCtr[parent][powLevel] ++;
            } //trapezoidal-wave
        } //pow
    } //level

}

//direction > 0 : triangle
//else : upside-down triangle
void FuncManager::powerCallHopelessReminder(int pivotLevel, int startPower, int endPower, int direction, int threadPerNode, const std::vector<int> *levelPtr, int numaLocalArg, int offset, int parent)
{
    int tid = omp_get_thread_num();
    int localTid = tid % threadPerNode;

    int curPow = startPower;

    if(direction > 0)
    {
        for(int base_len=(endPower-startPower); base_len>=0; --base_len)
        {
//Barrier has to be replaced with p2p, because only barrier within a node
#pragma omp barrier
            //reminder-head
            int startLevel = pivotLevel-base_len;
            int endLevel = pivotLevel+base_len;
            for(int level=startLevel; level<=endLevel; ++level)
            {
                SPLIT_LEVEL_PER_THREAD(level);
                //can be a function ptr
                powerFunc(startRow_tid, endRow_tid, curPow, numaLocalArg, args);
            }

            //TODO: make it general
            UNLOCK(endLevel-1, curPow+1);
            lockCtr[parent][endLevel] = threadPerNode*curPow;
            ++curPow;
        }
    }
    else
    {
        for(int base_len=0; base_len<=(startPower-endPower); ++base_len)
        {
#pragma omp barrier
            //reminder-head
            int startLevel = pivotLevel-base_len;
            int endLevel = pivotLevel+base_len;
            for(int level=startLevel; level<=endLevel; ++level)
            {
                SPLIT_LEVEL_PER_THREAD(level);
                //can be a function ptr
                powerFunc(startRow_tid, endRow_tid, curPow, numaLocalArg, args);
            }
            ++curPow;
        }
    }
}

void FuncManager::powerCallNodeReminder(int threadPerNode, int startSlope, int endSlope, const std::vector<int> *levelPtr, const std::vector<int> *nodePtr, int numaLocalArg, int offset)
{
    int tid = omp_get_thread_num();
    int localTid = tid % threadPerNode;
    int nodeGroup = tid / threadPerNode;

    for(int pow=1; pow<power; ++pow)
    {
#pragma omp barrier
        //reminder-head
        if(startSlope!=0)
        {
            int startLevel = nodePtr->at(nodeGroup);
            int endLevel = std::min(startLevel+pow, nodePtr->at(nodeGroup+1));
            for(int level=startLevel; level<endLevel; ++level)
            {
                SPLIT_LEVEL_PER_THREAD(level);
                //can be a function ptr
                powerFunc(startRow_tid, endRow_tid, pow+1, numaLocalArg, args);
            }
        }
        //reminder-tail
        if(endSlope!=0)
        {
            int endLevel = nodePtr->at(nodeGroup+1);
            int startLevel = std::max(nodePtr->at(nodeGroup), endLevel-pow);
            for(int level=startLevel; level<endLevel; ++level)
            {
                SPLIT_LEVEL_PER_THREAD(level);
                //can be a function ptr
                powerFunc(startRow_tid, endRow_tid, pow+1, numaLocalArg, args);
            }
        }
    }
}

void recursivePowerCall(FuncManager* funMan, int parent)
{
#pragma omp parallel
    {
        recursivePowerCallSerial(funMan, parent);
    }
}

void recursivePowerCallSerial(FuncManager* funMan, int parent)
{
//#pragma omp parallel
    {
        //printf("parent = %d\n", parent);
        funMan->initPowerRun(parent);
        //*nodePtr - for NUMA
        //**levelPtr - for power
        //int totalLevel = matPower->getTotalLevel();
        //do pointers to avoid copying expense; since tree is now stabe (size fixed)
        //this shouldn't be a problem
        std::vector<MPLeaf>* tree = &(funMan->matPower->tree);
        std::vector<int>* nodePtr = &(tree->at(parent).nodePtr);
        int totalNodes = nodePtr->size()-1;
        std::vector<int>* levelPtr  = &(tree->at(parent).levelPtr);
        std::vector<int>* unlockRow = &(tree->at(parent).unlockRow);
        std::vector<int>* dangerRow = &(tree->at(parent).dangerRow);
        std::vector<int>* unlockCtr = &(tree->at(parent).unlockCtr);
        std::vector<int>* unitPtr = &(tree->at(parent).unitPtr);
        std::vector<int>* unitNodePtr = &(tree->at(parent).unitNodePtr);

        //all delete after putting reminder with general
        //int power = funMan->power;
        bool numaSplit = funMan->numaSplit;
        int threadPerNode = funMan->activethreads/totalNodes;
        //printf("####### nthreads = %d\n", omp_get_num_threads());
        int tid = omp_get_thread_num();
        //int localTid = tid % threadPerNode;
        int nodeGroup = tid / threadPerNode;
        int numaLocalArg = (numaSplit)?nodeGroup:0;
        int offset = 0;
        offset = (parent==0)?levelPtr->at(nodePtr->at(numaLocalArg)):0;
        // printf("Offset = %d\n", offset);
        int startNode = unitNodePtr->at(nodeGroup);
        int endNode = unitNodePtr->at(nodeGroup+1);
        int unitCtr = startNode;
        int childCount = 0; //TODO: child count start for each node

        int nodeStartSlope=0, nodeEndSlope=0;
        //find start and end slope
        if(nodeGroup != 0)
        {
            nodeStartSlope = 1;
        }
        if(nodeGroup != (totalNodes-1))
        {
            nodeEndSlope = -1;
        }

        while(unitCtr < endNode)
        {
            int startSlope=-1, endSlope=-1;
            int startSkew = (funMan->power-1);
            //find start and end slope
            if(unitCtr == 0)
            {
                startSlope = nodeStartSlope;
                startSkew = 0;
            }
            if(unitCtr == (endNode-2))
            {
                endSlope = nodeEndSlope;
            }

            int startLevel = unitPtr->at(unitCtr)+startSkew;
            int endLevel = unitPtr->at(unitCtr+1);
            //printf("unitCtr = %d, startSlope = %d, endSlope = %d, startLevel = %d, endLevel = %d, startRow = %d, endRow = %d\n", unitCtr, startSlope, endSlope, startLevel, endLevel, levelPtr->at(startLevel), levelPtr->at(endLevel));
            //main-body
            funMan->powerCallGeneral(startLevel, endLevel, startLevel, endLevel, startSlope, endSlope, threadPerNode, levelPtr, unlockRow, unlockCtr, dangerRow, numaLocalArg, offset, parent);
            //printf("main finished\n");
            int hopelessEnd = unitPtr->at(unitCtr+2);
            if(hopelessEnd > endLevel)
            {
                //right-reminder of Hopeless, TODO currently only for power=2
                funMan->powerCallHopelessReminder(hopelessEnd, 1, 1, 1, threadPerNode, levelPtr, numaLocalArg, offset, parent);
                //printf("right-reminder finished\n");
                //do recursive call to hopeless region
                recursivePowerCallSerial(funMan, tree->at(parent).children[childCount]);
                //printf("recursive finished\n");
                ++childCount;
                //left-reminder of Hopeless, TODO currently only for power=2
                funMan->powerCallHopelessReminder(endLevel-1, 2, 2, -1, threadPerNode, levelPtr, numaLocalArg, offset, parent);
                //printf("left-reminder finished\n");
            }
            unitCtr = unitCtr+2;
        }

        //treat NUMA-interfaces
        if(totalNodes > 1)
        {
            funMan->powerCallNodeReminder(nodeStartSlope, nodeEndSlope, threadPerNode, levelPtr, nodePtr, numaLocalArg, offset);
        }

    } //parallel
}

#else

void FuncManager::powerRun()
{
    //int totalLevel = matPower->getTotalLevel();
    int totalNodes = matPower->getTotalNodes();
    int* levelPtr = matPower->getLevelPtrRef();
    int* nodePtr = matPower->getNodePtrRef();

#pragma omp parallel
    {
        int threadPerNode = omp_get_num_threads()/totalNodes;
        int tid = omp_get_thread_num();
        int nodeGroup = tid / threadPerNode;
        int localTid = tid % threadPerNode;
        //body
        {
            //printf("here is %d\n", sched_getcpu());
            int startLevel = nodePtr[nodeGroup];
            int endLevel = nodePtr[nodeGroup+1];
            int maxLevelCount = 0;
            for(int i=0; i<totalNodes; ++i)
            {
                maxLevelCount = std::max(maxLevelCount, (nodePtr[i+1]-nodePtr[i]));
            }
            int maxEndLevel = startLevel + maxLevelCount; //needed so that everyone calls barrier
            for(int level=startLevel; level<(maxEndLevel); ++level)
            {
                if(level < endLevel)
                {
                    for(int pow=0; pow<power; ++pow)
                    {
                        int powLevel = (level-pow);

                        if( (powLevel >= (startLevel+pow)) && (powLevel < (endLevel-pow)) )
                        {
                            SPLIT_LEVEL_PER_THREAD(powLevel);
                            //can be a function ptr
                            powerFunc(startRow_tid, endRow_tid, pow+1, args);
#if 0
#pragma omp parallel for schedule(static)
                            for(int row=levelPtr[powLevel]; row<levelPtr[powLevel+1]; ++row)
                            {
                                double tmp = 0;
                                for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
                                {
                                    tmp += A[idx]*x[(pow)*graph->NROW+col[idx]];
                                }
                                x[(pow+1)*graph->NROW+row] = tmp;
                            }
#endif
                        }
#pragma omp barrier
                        //for correctness it has to be here, but its rarely that this conflict happens
                        //to detect this run for example FDM-512 with p=80
                        //I think we need this actually, but this actually
                        //separates into inCache and burst phase
                    }
                }
                //#pragma omp barrier

            }
        }
        //#pragma omp barrier
        //reminder
        for(int pow=1; pow<power; ++pow)
        {
            //reminder-head
            {
                int startLevel = nodePtr[nodeGroup];
                int endLevel = std::min(startLevel+pow, nodePtr[nodeGroup+1]);
                for(int level=startLevel; level<endLevel; ++level)
                {
                    SPLIT_LEVEL_PER_THREAD(level);
                    //can be a function ptr
                    powerFunc(startRow_tid, endRow_tid, pow+1, args);
                }
            }
            //reminder-tail
            {
                int endLevel = nodePtr[nodeGroup+1];
                int startLevel = std::max(nodePtr[nodeGroup], endLevel-pow);
                for(int level=startLevel; level<endLevel; ++level)
                {
                    SPLIT_LEVEL_PER_THREAD(level);
                    //can be a function ptr
                    powerFunc(startRow_tid, endRow_tid, pow+1, args);
                }

            }
#pragma omp barrier
        }
    }
}

#endif

void FuncManager::Run(bool rev_)
{
    if(!power_fn)
    {
        rev = rev_;

        if(zoneTree == NULL)
        {
            ERROR_PRINT("NO Zone Tree present; Have you registered the function");
        }

        /*        if(rev)
                  {
                  func(serialPart[1],serialPart[0],args);
                  }
                  */
        int root = 0;
#ifdef RACE_KERNEL_THREAD_OMP
        int resetNestedState = omp_get_nested();
        int resetDynamicState = omp_get_dynamic();
        //set nested parallelism
        omp_set_nested(1);
        omp_set_dynamic(0);
        recursiveFun(root);

        //reset states
        omp_set_nested(resetNestedState);
        omp_set_dynamic(resetDynamicState);
#else

        // START_TIME(main_fn_call);
        /*recursiveCall(this, root);*/
        recursiveFun(root);
        //STOP_TIME(main_fn_call)
#endif

        if(!rev)
        {
            func(serialPart[0],serialPart[1],args);
        }
        if(rev)
        {
            func(serialPart[1],serialPart[0],args);
        }


    }
    else
    {
#ifdef RACE_KERNEL_THREAD_OMP
        int resetNestedState = omp_get_nested();
        //set nested parallelism
        //printf("setting nested\n");
        omp_set_nested(0);
        recursivePowerFun(0);
        //powerRun();

        //reset states
        omp_set_nested(resetNestedState);
#endif

    }
}

/*void FuncManager::RunOMP()
  {
//test_omp(a,b,c,d,0,len,1);
}*/

bool FuncManager::isNumaSplit()
{
    return numaSplit;
}
