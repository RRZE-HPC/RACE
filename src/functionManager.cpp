#include "functionManager.h"
#include "test.h"

FuncManager::FuncManager(void (*f_) (int,int,void*), void*args_, ZoneTree *zoneTree_, LevelPool* pool_):func(f_),args(args_),zoneTree(zoneTree_), pool(pool_)
{
//    func = std::bind(f_,*this,std::placeholders::_1,std::placeholders::_2,args);
/*    len = 1024*1024;
    a = new double[len];
    b = new double[len];
    c = new double[len];
    d = new double[len];*/
}

FuncManager::~FuncManager()
{
/*    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    */
}


#ifdef NAME_KERNEL_THREAD_OMP
//OMP version nested pinning not working
void FuncManager::recursiveCall(int parent)
{
    std::vector<int> *children = &(zoneTree->at(parent).childrenZ);
    int currNthreads = static_cast<int>(children->size()/2.0);
    if(currNthreads <= 1)
    {
        std::vector<int> range = zoneTree->at(parent).valueZ;
        func(range[0],range[1],args);
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
            recursiveCall(children->at(2*tid));
#pragma omp barrier
            //printf("Black: tid = %d, pinOrder = %d, cpu = %d\n",tid,pinOrder,sched_getcpu());
            //do Black
            recursiveCall(children->at(2*tid+1));
        }
    }
}
#else
//PTHREAD implementation
void FuncManager::recursiveCall(int parentIdx)
{
    std::vector<int>* children = &(zoneTree->at(parentIdx).childrenZ);
    int nthreads = children->size()/2;
    if(nthreads <= 1)
    {
        //        int pinOrder = zoneTree->at(parentIdx).pinOrder;
        //        printf(" pinOrder = %d, cpu = %d\n",pinOrder,sched_getcpu());
        std::vector<int>* range = &zoneTree->at(parentIdx).valueZ;
        func(range->at(0),range->at(1),args);
        //  test(1);
        // sleep(1);
    }
    else
    {
        for(int tid=0; tid<nthreads; ++tid)
        {
            //pool->tree[parentIdx].addJob(tid, std::bind(test,a,b,c,d, tid*len/2.0, (tid+1)*len/2.0, std::placeholders::_1), 1);
            pool->tree[parentIdx].addJob(tid, std::bind(&FuncManager::recursiveCall,this, std::placeholders::_1), children->at(2*tid) );
        }
        pool->tree[parentIdx].barrier();
        for(int tid=0; tid<nthreads; ++tid)
        {
            //pool->tree[parentIdx].addJob(tid, std::bind(test,a,b,c,d, tid*len/2.0, (tid+1)*len/2.0, std::placeholders::_1), 1);
            pool->tree[parentIdx].addJob(tid, std::bind(&FuncManager::recursiveCall,this, std::placeholders::_1), children->at(2*tid+1));
        }
        pool->tree[parentIdx].barrier();
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

    recursiveCall(root);

    //reset states
    omp_set_nested(resetNestedState);
    omp_set_dynamic(resetDynamicState);
#else
    START_TIME(main_fn_call);
    recursiveCall(root);
    STOP_TIME(main_fn_call)
#endif
}

