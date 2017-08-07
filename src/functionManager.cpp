#include "functionManager.h"
#include "test.h"
#include <iostream>
#include <typeinfo>

FuncManager::FuncManager(void (*f_) (int,int,void*), void*args_, ZoneTree *zoneTree_, LevelPool* pool_):func(f_),args(args_),zoneTree(zoneTree_), pool(pool_)
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
    printf("Copy constructor called\n");
    func = obj.func;
    std::cout<<"fn type = "<<typeid(obj.func).name() << std::endl;
    printf("copied fn\n");
    zoneTree = obj.zoneTree;
    args = obj.args;
    pool = obj.pool;
    printf("finished copying\n");
}

FuncManager::~FuncManager()
{
    /*  delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
        */
}


#ifdef NAME_KERNEL_THREAD_OMP
//OMP version nested pinning not working
void recursiveCall(FuncManager* funMan, int parent)
{
    std::vector<int> *children = &(funMan->zoneTree->at(parent).childrenZ);
    int blockPerThread = getBlockPerThread(funMan->zoneTree->dist, funMan->zoneTree->d2Type);
    int currNthreads = static_cast<int>(children->size()/blockPerThread);
    if(currNthreads <= 1)
    {
        std::vector<int> range = funMan->zoneTree->at(parent).valueZ;
        funMan->func(range[0],range[1],funMan->args);
    }
    else
    {
#pragma omp parallel num_threads(currNthreads)
        {
            int tid = omp_get_thread_num();
            //Pin in each call
            //int pinOrder = zoneTree->at(children->at(2*tid)).pinOrder;
            //printf("omp_proc_bind = %d\n", omp_get_proc_bind());
            //pool->pin.pinThread(pinOrder);
            //printf("child = %d, Red: tid = %d, pinOrder = %d, cpu = %d\n",children->at(2*tid),tid,pinOrder,sched_getcpu());
            //do Red
            for(int block=0; block<blockPerThread-1; ++block)
            {
                recursiveCall(funMan, children->at(blockPerThread*tid+block));
                #pragma omp barrier
            }
            //printf("Black: tid = %d, pinOrder = %d, cpu = %d\n",tid,pinOrder,sched_getcpu());
            //do Black
            recursiveCall(funMan, children->at(blockPerThread*(tid+1)-1));
        }
    }
}
#else
//PTHREAD implementation
void recursiveCall(FuncManager* funMan, int parentIdx)
{
    std::vector<int>* children = &(funMan->zoneTree->at(parentIdx).childrenZ);
    int blockPerThread = getBlockPerThread(funMan->zoneTree->dist, funMan->zoneTree->d2Type);
    int nthreads = children->size()/blockPerThread;
    if(nthreads <= 1)
    {
        //        int pinOrder = zoneTree->at(parentIdx).pinOrder;
        //        printf(" pinOrder = %d, cpu = %d\n",pinOrder,sched_getcpu());
        std::vector<int>* range = &funMan->zoneTree->at(parentIdx).valueZ;
        START_TIME(func);
        funMan->func(range->at(0),range->at(1), funMan->args);
        STOP_TIME(func);
        funMan->zoneTree->at(parentIdx).time += GET_TIME(func);
        //test(funMan->a,funMan->b,funMan->c,funMan->d,range->at(0), range->at(1),1);
        //  test(1);
        // sleep(1);
    }
    else
    {
        for(int block=0; block<blockPerThread; ++block)
        {
            for(int tid=0; tid<nthreads; ++tid)
            {
                //pool->tree[parentIdx].addJob(tid, std::bind(test,a,b,c,d, tid*len/2.0, (tid+1)*len/2.0, std::placeholders::_1), 1);
                funMan->pool->tree[parentIdx].addJob(tid, funMan->recursiveFun, children->at(blockPerThread*tid+block) );
            }
            funMan->pool->tree[parentIdx].barrier();
            /*if(parentIdx == 0)
              {
              funMan->barrierTime=funMan->pool->tree[parentIdx].barrierTime;
              }*/
        }
    }
}

#endif

void FuncManager::Run()
{
    if(zoneTree == NULL)
    {
        ERROR_PRINT("NO Zone Tree present; Have you registered the function");
    }

    int root = 0;
#ifdef NAME_KERNEL_THREAD_OMP
    int resetNestedState = omp_get_nested();
    int resetDynamicState = omp_get_dynamic();
    //set nested parallelism
    omp_set_nested(1);
    omp_set_dynamic(0);
    recursiveCall(this, root);

    //reset states
    omp_set_nested(resetNestedState);
    omp_set_dynamic(resetDynamicState);
#else

    // START_TIME(main_fn_call);
    /*recursiveCall(this, root);*/
    recursiveFun(root);
    //STOP_TIME(main_fn_call)
#endif
}

void FuncManager::RunOMP()
{
    //test_omp(a,b,c,d,0,len,1);
}

