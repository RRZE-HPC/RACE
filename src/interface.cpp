#include "interface.h"
#include "simdify.h"
#include <algorithm>
#include <iostream>

RACE::Interface::Interface(int nrow_,int nthreads_, RACE::dist dist_, int *rowPtr_, int *col_, int SMT_, RACE::PinMethod pinMethod_, int *initPerm_, int *initInvPerm_, RACE::d2Method d2Type_):graph(NULL),nrow(nrow_),distance(dist_),d2Type(d2Type_),requestedThreads(nthreads_),availableThreads(-1),SMT(SMT_),pinMethod(pinMethod_),pool(NULL),initPerm(initPerm_),initInvPerm(initInvPerm_),rowPtr(rowPtr_),col(col_),zoneTree(NULL)
{

}

RACE::Interface::~Interface()
{
    if(zoneTree) {
        delete zoneTree;
    }

    if(graph) {
        delete graph;
    }

    if(pool) {
        delete pool;
    }

    for(unsigned i=0; i<funMan.size(); ++i)
    {
        delete funMan[i];
    }
}


void RACE::Interface::RACEColor()
{

    //1. Construct Graph
    graph = new Graph(nrow, nrow, rowPtr, col, initPerm, initInvPerm);

    //2. Call level_recursion
    LevelRecursion lr(graph, requestedThreads, distance, d2Type);
    lr.levelBalancing();
    availableThreads = lr.getAvailableThreads();
    int len;
    lr.getPerm(&perm, &len);
    lr.getInvPerm(&invPerm, &len);
    zoneTree = lr.getZoneTree();

    pool = new LevelPool(zoneTree, SMT, pinMethod);

#ifdef RACE_KERNEL_THREAD_OMP
    pool->pin.pinApplication();
#else
    pool->createPool();//creates pinned thread pools
#endif
    printZoneTree();

    /*    printf("Checking Coloring\n");
          if(D2Checker())
          {
          ERROR_PRINT("Conflict in coloring\n");
          }
          printf("Checking Finished\n");
          */    
}


void RACE::Interface::printZoneTree()
{
   zoneTree->printTree();
}

void RACE::Interface::getPerm(int **perm_, int *len_)
{
    if(initPerm)
    {
        int *totPerm = new int [nrow];
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<nrow; ++i)
        {
            totPerm[i] = initPerm[perm[i]];
        }

        (*perm_) = totPerm;
        delete[] perm;
    }
    else
    {
        (*perm_) = perm;
    }
    (*len_) = nrow;
}

void RACE::Interface::getInvPerm(int **invPerm_, int *len_)
{
    if(initInvPerm)
    {
        int *totInvPerm = new int [nrow];
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<nrow; ++i)
        {
            totInvPerm[i] = invPerm[initInvPerm[i]];
        }


        (*invPerm_) = totInvPerm;
        delete[] invPerm;
    }
    else
    {
        (*invPerm_) = invPerm;
    }
    (*len_) = nrow;
}

int RACE::Interface::getNumThreads()
{
    return availableThreads;
}


int RACE::Interface::registerFunction(void (*f) (int,int,void *), void* args)
{
    //    pool->pin.pinApplication();

    FuncManager* currFun = new FuncManager(f, args, zoneTree, pool);
    funMan.push_back(currFun);

    //     funMan = new FuncManager(f, args, zoneTree, pool);

    /* for(int i=0; i<10; ++i) {
       START_TIME(omp_fn);
       funMan->RunOMP();
       STOP_TIME(omp_fn);
       }*/
    //TODO just pin for the first time; but now OpenMP does not work with this
    //model
    //Pin threads
    //Pin pin(zoneTree, true);
    //pin.pinApplication();

    return (funMan.size()-1);
}

void RACE::Interface::executeFunction(int funcId)
{
    funMan[funcId]->Run();
    //  funMan->Run();
}


void RACE::Interface::resetTime()
{
    zoneTree->resetTime();
}


bool RACE::Interface::detectConflict(std::vector<int> range1, std::vector<int> range2)
{
    bool conflict = false;
    for(int i=range1[0]; i<range1[1]; ++i)
    {
        std::vector<int>* children_i = &(graph->at(i).children);
        for(int j=range2[0]; j<range2[1]; ++j)
        {
            std::vector<int>* children_j = &(graph->at(j).children);
            conflict = (std::find_first_of(children_i->begin(), children_i->end(), children_j->begin(), children_j->end()) != children_i->end());
            if(conflict)
            {
                ERROR_PRINT("Conflict at row: %d %d\nCheck i: [",i,j);
                for(unsigned child_i=0; child_i<children_i->size(); ++child_i)
                {
                    printf("%d, ",children_i->at(child_i));
                }
                printf("]\nCheck j: [");
                for(unsigned child_j=0; child_j<children_j->size(); ++child_j)
                {
                    printf("%d, ",children_j->at(child_j));
                }
                printf("]\n");

                auto result = std::find_first_of(children_i->begin(), children_i->end(), children_j->begin(), children_j->end());
                std::cout<<"conflicted values = " << std::distance(children_i->begin(),result)<< "and "<< std::distance(children_j->begin(),result)<<std::endl;

                break;
            }
        }
        if(conflict)
        {
            break;
        }
    }

    return conflict;
}

bool RACE::Interface::recursiveChecker(int parent)
{
    std::vector<int> *children = &(zoneTree->at(parent).children);
    int totalSubBlocks = zoneTree->at(parent).totalSubBlocks;
    int blockPerThread = getBlockPerThread(distance, d2Type);
    bool conflict = false;

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        int currNThreads = (children->at(2*subBlock+1)-children->at(2*subBlock))/blockPerThread;
        printf("checking idx = %d\n",parent);
        if(currNThreads > 1)
        {
            for(int start=0; start<blockPerThread; ++start)
            {
                for(int i=start; i<currNThreads; i+=blockPerThread)
                {
                    std::vector<int> rangeOuter = zoneTree->at(children->at(2*subBlock)+i).valueZ;
                    for(int j=start; j<=currNThreads; j+=2)
                    {
                        if(i!=j)
                        {
                            std::vector<int> rangeInner = zoneTree->at(children->at(2*subBlock)+j).valueZ;
                            conflict = detectConflict(rangeOuter, rangeInner);
                            if(conflict)
                            {
                                break;
                            }
                        }
                    }
                    if(conflict)
                    {
                        break;
                    }
                }
            }
            for(int child=0; child< blockPerThread*currNThreads; ++child)
            {
                recursiveChecker(children->at(2*subBlock)+child);
            }
        }
    }
    return conflict;
}

bool RACE::Interface::D2Checker()
{
    bool conflict = false;
    conflict = recursiveChecker(5);

    return conflict;
}

void RACE::Interface::sleep()
{
    pool->sleepPool();
}

bool RACE::Interface::simdify(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* clp, double* val, bool diagFirst)
{
    return simdifyTemplate<double> (simdWidth, C, nrows, col_new, chunkStart, rl, clp, val, this, diagFirst);
}

bool RACE::Interface::simdify(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* clp, float* val, bool diagFirst)
{
    return simdifyTemplate<float> (simdWidth, C, nrows, col_new, chunkStart, rl, clp, val, this, diagFirst);
}

bool RACE::Interface::simdifyD1(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* cl, double* val, bool d2Compatible)
{
    return simdifyD1Template<double> (simdWidth, C, nrows, col_new, chunkStart, rl, cl, val, d2Compatible);
}

bool RACE::Interface::simdifyD1(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* cl, float* val, bool d2Compatible)
{
    return simdifyD1Template<float> (simdWidth, C, nrows, col_new, chunkStart, rl, cl, val, d2Compatible);
}


void RACE::Interface::pinThread(int threadId)
{
    pool->pin.pinThread(threadId);
}
