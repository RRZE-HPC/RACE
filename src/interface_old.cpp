#include "interface.h"

RACEInterface::RACEInterface(int nrow_,int nthreads_, RACE::dist dist_, int *rowPtr_, int *col_, int *initPerm_, int *initInvPerm_):nrow(nrow_),dist(dist_),requestedThreads(nthreads_),availableThreads(-1),initPerm(initPerm_),initInvPerm(initInvPerm_),rowPtr(rowPtr_),col(col_)
{
}

void RACEInterface::RACEColor()
{
    //1. Construct Graph
    Graph bmc(nrow, nrow, rowPtr, col, initPerm, initInvPerm);

    //2. Traversal
    Traverse traverser(&bmc, TWO);
    traverser.calculateDistance();
    LevelData *levelData;
    levelData = traverser.getLevelData();
    traverser.getPerm(&perm, &permLen);
    traverser.getInvPerm(&invPerm, &invPermLen);

    //3. Load Balncing
    LB lb(requestedThreads, levelData);
    lb.D2LB();
    lb.getZonePtr(&zonePtr, &zonePtrLen);
    availableThreads = lb.getNumThreads();
    delete levelData;
}

void RACEInterface::getZonePtr(int **zonePtr_, int *len_)
{
    (*zonePtr_) = zonePtr;
    (*len_) = zonePtrLen;
}

void RACEInterface::getPerm(int **perm_, int *len_)
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
    (*len_) = permLen;
}

void RACEInterface::getInvPerm(int **invPerm_, int *len_)
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
    (*len_) = invPermLen;
}

int RACEInterface::getNumThreads()
{
    return availableThreads;
}


