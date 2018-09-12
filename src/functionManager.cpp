#include "functionManager.h"
#include "test.h"
#include <iostream>
#include <typeinfo>

FuncManager::FuncManager(void (*f_) (int,int,void*), void *args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart_):rev(false),func(f_),args(args_),zoneTree(zoneTree_), pool(pool_), serialPart(serialPart_)
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
}

FuncManager::FuncManager(const FuncManager &obj)
{
    rev = obj.rev;
    func = obj.func;
    zoneTree = obj.zoneTree;
    args = obj.args;
    pool = obj.pool;
    serialPart = obj.serialPart;
}

FuncManager::~FuncManager()
{
    /*  delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
        */
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



void FuncManager::Run(bool rev_)
{
    rev = rev_;

    if(zoneTree == NULL)
    {
        ERROR_PRINT("NO Zone Tree present; Have you registered the function");
    }

    if(rev)
    {
        func(serialPart[1],serialPart[0],args);
    }

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
}

/*void FuncManager::RunOMP()
{
    //test_omp(a,b,c,d,0,len,1);
}*/

