#include "functionManager.h"

FuncManager::FuncManager(void (*f_) (int,int,void*), void*args_, ZoneTree *zoneTree_, Pin* pin_):func(f_),args(args_),zoneTree(zoneTree_), pin(pin_)
{
//    func = std::bind(f_,*this,std::placeholders::_1,std::placeholders::_2,args);
}

FuncManager::~FuncManager()
{

}

//OMP nested pinning not working
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
            int pinOrder = zoneTree->at(children->at(2*tid)).pinOrder;
            //Pin in each call
            pin->pinThread(pinOrder);
            //printf("Red: tid = %d, pinOrder = %d, cpu = %d\n",tid,pinOrder,sched_getcpu());
            //do Red
            recursiveCall(children->at(2*tid));
#pragma omp barrier
            //printf("Black: tid = %d, pinOrder = %d, cpu = %d\n",tid,pinOrder,sched_getcpu());
            //do Black
            recursiveCall(children->at(2*tid+1));
        }
    }
}


void FuncManager::Run()
{
    if(zoneTree == NULL)
    {
        ERROR_PRINT("NO Zone Tree present; Have you registered the function");
    }

//#ifdef NAME_HAVE_OPENMP
    int resetNestedState = omp_get_nested();
    int resetDynamicState = omp_get_dynamic();
    //set nested parallelism
    omp_set_nested(1);
    omp_set_dynamic(0);

    int root = 0;
    recursiveCall(root);

    //reset states
    omp_set_nested(resetNestedState);
    omp_set_dynamic(resetDynamicState);
//#endif
}

