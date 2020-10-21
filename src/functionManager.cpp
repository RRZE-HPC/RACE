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
FuncManager::FuncManager(void (*f_) (int,int,void*), void *args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart_):rev(false),power_fn(false), numaSplit(false), func(f_),args(args_),zoneTree(zoneTree_), pool(pool_), serialPart(serialPart_), lockCtr(NULL)
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
FuncManager::FuncManager(void (*f_) (int,int,int,int,void*), void *args_, int power_, mtxPowerRecursive *matPower_, bool numaSplit_):power_fn(true), numaSplit(numaSplit_), powerFunc(f_), args(args_), power(power_), matPower(matPower_), lockCtr(NULL), lockTableCtr(NULL)
{
#ifdef POWER_WITH_FLUSH_LOCK
    initPowerRun();
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
        if(lockCtr)
        {
            delete[] lockCtr;
        }
#ifdef POWER_WITH_FLUSH_LOCK
        else
        {
            ERROR_PRINT("Lock counter does not exist");
        }
#endif
    }
    if(power_fn)
    {
        if(lockTableCtr)
        {
            free((void*)lockTableCtr);
        }
#ifdef POWER_WITH_FLUSH_LOCK
        else
        {
            ERROR_PRINT("Lock table counter does not exist");
        }
#endif
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

void FuncManager::initPowerRun()
{
    //need to move this routine to matrixPower.cpp
#ifdef POWER_WITH_FLUSH_LOCK
    std::vector<int> nodePtr = matPower->tree[0].nodePtr;
    int totalNodes = nodePtr.size()-1;
    int totalLevels = matPower->tree[0].levelPtr.size()-1;

    //allocate only if needed
    if(!lockCtr)
    {
        lockCtr = new volatile int[totalLevels];
    }
#pragma omp parallel
    {
        int threadPerNode = omp_get_num_threads()/totalNodes;
        int tid = omp_get_thread_num();
        int nodeGroup = tid / threadPerNode;
        int localTid = tid % threadPerNode;
        //printf("here is %d\n", sched_getcpu());
        int startLevel = nodePtr[nodeGroup];
        int endLevel = nodePtr[nodeGroup+1];

        //only 1 thread per group, maintain NUMA
        if(localTid == 0)
        {
            for(int level=startLevel; level<endLevel; ++level)
            {
                lockCtr[level] = 0;
            }
        }
    }

    //pow-level semi locks
    if(lockTableCtr)
    {
        free((void*)lockTableCtr);
    }
    lockTableCtr = (volatile int*)malloc(sizeof(volatile int)*power*totalLevels);

    for(int l=0; l<totalLevels; ++l)
    {
        for(int p=0; p<power; ++p)
        {
            lockTableCtr[l*power+p] = 0;
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
        lockTableCtr[_level_*power+_pow_]++;\
    }\
}

#define WAIT(_level_, _pow_)\
{\
    if(_pow_ != 0)\
    {\
        while(lockTableCtr[_level_*power+_pow_] != unlockCtr[_level_+1])\
        {\
            _mm_pause();\
            /*printf("Waiting on %d, need to reach%d, pow = %d, level = %d, startRow_tid = %d, unlock = %d, danger = %d, endRow_tid = %d\n", lockTableCtr[_level_*power+_pow_], unlockCtr[_level_], _pow_, _level_, startRow_tid, currUnlockRow, dangerRowStart, endRow_tid);*/\
        }\
    }\
}\


void recursivePowerCall(FuncManager* funMan, int parent)
{

    funMan->initPowerRun();
    //*nodePtr - for NUMA
    //**levelPtr - for power
    //int totalLevel = matPower->getTotalLevel();
    //do pointers to avoid copying expense; since tree is now stabe (size fixed)
    //this shouldn't be a problem
    std::vector<MPLeaf>* tree = &(funMan->matPower->tree);
    std::vector<int> nodePtr = tree->at(parent).nodePtr;
    int totalNodes = nodePtr.size()-1;
    std::vector<int> levelPtr = tree->at(parent).levelPtr;
    std::vector<int> unlockRow = tree->at(parent).unlockRow;
    std::vector<int> dangerRow = tree->at(parent).dangerRow;
    std::vector<int> unlockCtr = tree->at(parent).unlockCtr;

    int power = funMan->power;
    int activethreads = funMan->activethreads;
    bool numaSplit = funMan->numaSplit;
    void *args = funMan->args;
    volatile int* lockCtr = funMan->lockCtr;
    volatile int* lockTableCtr = funMan->lockTableCtr;
    //don't spawn threads recursively; TODO switch off recursive before calling this function
#pragma omp parallel
    {
        int threadPerNode = activethreads/totalNodes;
        int tid = omp_get_thread_num();
        int nodeGroup = tid / threadPerNode;
        int localTid = tid % threadPerNode;
        int numaLocalArg = (numaSplit)?nodeGroup:0;
        int offset = levelPtr[nodePtr[numaLocalArg]];
        // printf("Offset = %d\n", offset);
        int startLevel = nodePtr[nodeGroup];
        int endLevel = nodePtr[nodeGroup+1];

        //body
        {
            //printf("here is %d\n", sched_getcpu());
            /* int maxLevelCount = 0;
               for(int i=0; i<totalNodes; ++i)
               {
               maxLevelCount = std::max(maxLevelCount, (nodePtr->at(i+1)-nodePtr->at(i)));
               }
               int maxEndLevel = startLevel + maxLevelCount; //needed so that everyone calls barrier
               */
            for(int level=startLevel; level<(endLevel); ++level)
            {
                //tid responsible for (level-1) lock, locks the level
                for(int pow=0; pow<power; ++pow)
                {
                    int powLevel = (level-pow);

                    //square wave
                    if( ((powLevel >= (startLevel)) && (powLevel < (endLevel))) )
                    {
                        //trapezoidal wave
                        if( ((powLevel >= (startLevel+pow)) && (powLevel < (endLevel-pow))) )
                        {

                            //wait till powLevel is free; use ctr[powLevel] == threads*(pow)
                            while(lockCtr[powLevel] < threadPerNode*pow)
                            {
                                _mm_pause();
#ifdef RACE_DEBUG
                                printf("lock ctr = %d, pow = %d, level = %d\n", lockCtr[powLevel], pow, powLevel);
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
                                funMan->powerFunc(startRow_tid, endRow_tid, pow+1, numaLocalArg, args);
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
                                    funMan->powerFunc(startRow_tid, till_row, pow+1, numaLocalArg, args);
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
                                //funMan->powerFunc(startRow_tid, endRow_tid, pow+1, args);
                                //till_row = endRow_tid;
                                }*/
                                if(dangerRowStart > endRow_tid)
                                {
                                    funMan->powerFunc(till_row, endRow_tid, pow+1, numaLocalArg, args);
#ifdef RACE_DEBUG
                                    printf("2. call tid = %d, [%d, %d], pow = %d, level = %d\n", omp_get_thread_num(), till_row, endRow_tid, pow, powLevel);
#endif
                                }
                                else
                                {
                                    funMan->powerFunc(till_row, dangerRowStart, pow+1, numaLocalArg, args);
#ifdef RACE_DEBUG
                                    printf("3. call tid = %d, [%d, %d], pow = %d, level = %d\n", omp_get_thread_num(), till_row, dangerRowStart, pow, powLevel);
#endif

                                    //check lock
                                    WAIT(powLevel, pow);

                                    funMan->powerFunc(dangerRowStart, endRow_tid, pow+1, numaLocalArg, args);
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
                            lockCtr[powLevel] ++;
                        }
                        else //trapezoidal-wave; only unlocking is done, no computation, computation done in reminder loop
                        {


                            SPLIT_LEVEL_PER_THREAD_P2P(powLevel);
#ifdef RACE_DEBUG
                            printf("** power = %d, pow = %d, level = %d, tid = %d, start row = %d, unlock row = %d, danger row = %d, end_row = %d\n", power, pow, powLevel, omp_get_thread_num(), startRow_tid, currUnlockRow, dangerRowStart, endRow_tid);
#endif
                            if(0)
                            {

                                ERROR_PRINT("Should never be here");
                                //To suppress unused variable warning
                                printf("** power = %d, pow = %d, level = %d, tid = %d, start row = %d, unlock row = %d, danger row = %d, end_row = %d\n", power, pow, powLevel, omp_get_thread_num(), startRow_tid, currUnlockRow, dangerRowStart, endRow_tid);
                            }
                            //if( (startRow_tid < currUnlockRow) && (endRow_tid >= currUnlockRow) ) //only 1 thread will be here
                            //                             #pragma omp single
                            if(startRow_tid < currUnlockRow) //only 1 thread will be here
                            {
                                //unlock powLevel -1 for pow+1
                                //int z = (pow) + (powLevel);
                                //int cons_term = std::max(0,z-power); //the term for having consequtive numbers
                                //int powLevel_to_unlock = ((z*(z+1))>>1) + (pow+1) - (cons_term*(cons_term+1)>>1);

                                UNLOCK((powLevel-1), (pow+1));
#ifdef RACE_DEBUG
                                printf("tid = %d, Unlocking Level = %d, pow = %d\n", omp_get_thread_num(), powLevel-1, pow+1);
#endif
                            }

                        } //trapezoidal-wave
                    }//square-wave
                } //pow
            } //level
        }//body
        //reminder
        for(int pow=1; pow<power; ++pow)
        {
#pragma omp barrier
            //reminder-head
            {
                int startLevel = nodePtr[nodeGroup];
                int endLevel = std::min(startLevel+pow, nodePtr[nodeGroup+1]);
                for(int level=startLevel; level<endLevel; ++level)
                {
                    SPLIT_LEVEL_PER_THREAD(level);
                    //can be a function ptr
                    funMan->powerFunc(startRow_tid, endRow_tid, pow+1, numaLocalArg, args);
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
                    funMan->powerFunc(startRow_tid, endRow_tid, pow+1, numaLocalArg, args);
                }

            }
        } //pow
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
