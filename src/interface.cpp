#include "interface.h"
#include <algorithm>
#include <iostream>

NAMEInterface::NAMEInterface(int nrow_,int nthreads_, dist_t dist_, int *rowPtr_, int *col_, int SMT_, PinMethod pinMethod_, int *initPerm_, int *initInvPerm_):graph(NULL),nrow(nrow_),dist(dist_),requestedThreads(nthreads_),availableThreads(-1),SMT(SMT_),pinMethod(pinMethod_),pool(NULL),initPerm(initPerm_),initInvPerm(initInvPerm_),rowPtr(rowPtr_),col(col_),zoneTree(NULL)
{

}

NAMEInterface::~NAMEInterface()
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


void NAMEInterface::NAMEColor()
{

    //1. Construct Graph
    graph = new Graph(nrow, nrow, rowPtr, col, initPerm, initInvPerm);

    //2. Call level_recursion
    LevelRecursion lr(graph, requestedThreads, dist, TWO_BLOCK);
    lr.levelBalancing();
    availableThreads = lr.getAvailableThreads();
    int len;
    lr.getPerm(&perm, &len);
    lr.getInvPerm(&invPerm, &len);
    zoneTree = lr.getZoneTree();

    pool = new LevelPool(zoneTree, SMT, pinMethod);

#ifdef NAME_KERNEL_THREAD_OMP
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


void NAMEInterface::printZoneTree()
{
    zoneTree->updateTime();
    for(unsigned i=0; i<zoneTree->tree->size(); ++i)
    {
        ZoneLeaf* currLeaf = &(zoneTree->at(i));
        printf("%d. Range:[", i);

        for(unsigned j=0; j<currLeaf->valueZ.size(); ++j)
        {
            printf("%d, ",currLeaf->valueZ[j]);
        }

        printf("] Children:[");


        for(unsigned j=0; j<currLeaf->childrenZ.size(); ++j)
        {
            printf("%d, ",currLeaf->childrenZ[j]);
        }
        printf("], Parent: %d, nthreads: %d, idealNThreads:%d, effRow: %d reachedLimit: %s, pinOrder = %d funTime = %f\n", currLeaf->parentZ, currLeaf->nthreadsZ, currLeaf->idealNthreadsZ, currLeaf->effRowZ, (currLeaf->reachedLimit)?"true":"false", currLeaf->pinOrder, currLeaf->time);
    }
}

void NAMEInterface::getPerm(int **perm_, int *len_)
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

void NAMEInterface::getInvPerm(int **invPerm_, int *len_)
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

int NAMEInterface::getNumThreads()
{
    return availableThreads;
}


int NAMEInterface::registerFunction(void (*f) (int,int,void *), void* args)
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

void NAMEInterface::executeFunction(int funcId)
{
    funMan[funcId]->Run();
    //  funMan->Run();
}


void NAMEInterface::resetTime()
{
    zoneTree->resetTime();
}


bool NAMEInterface::detectConflict(std::vector<int> range1, std::vector<int> range2)
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

bool NAMEInterface::recursiveChecker(int parent)
{
    std::vector<int> *children = &(zoneTree->at(parent).childrenZ);
    int currNThreads = static_cast<int>(children->size()/2.0);
    bool conflict = false;
    printf("checking idx = %d\n",parent);
    if(currNThreads > 1)
    {
        for(int start=0; start<2; ++start)
        {
            for(int i=start; i<=currNThreads; i+=2)
            {
                std::vector<int> rangeOuter = zoneTree->at(children->at(i)).valueZ;
                for(int j=start; j<=currNThreads; j+=2)
                {
                    if(i!=j)
                    {
                        std::vector<int> rangeInner = zoneTree->at(children->at(j)).valueZ;
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
        for(unsigned child=0; child< children->size(); ++child)
        {
            recursiveChecker(children->at(child));
        }
    }
    return conflict;
}

bool NAMEInterface::D2Checker()
{
    bool conflict = false;
    conflict = recursiveChecker(5);

    return conflict;
}
