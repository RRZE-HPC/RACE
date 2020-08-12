#ifndef _MATRIX_POWER_RECURSIVE_H
#define _MATRIX_POWER_RECURSIVE_H

#include "matrixPower.h"

class mtxPowerRecursive
{
    Graph* graph;

    //final values
    int* levelPtr;
    int* nodePtr;
    int* unlockRow;//for p2p sync
    int* dangerRow;// for p2p sync
    int* unlockCtr;

    int totalRows;
    int totalLevel;
    int highestPower;
    int numSharedCache;
    double cacheSize;
    double safetyFactor;

        int* perm;
    int* invPerm;

    std::vector<int> cache_violation_cutoff;
    int get_cache_violation_cutoff(int stage);


    public:
    mtxPowerRecursive(Graph* graph_, int highestPower_, int numSharedCache, double cacheSize_, double safetyFactor_);
    ~mtxPowerRecursive();

    void findPartition();
    //these are final values
    void getPerm(int **perm, int *len);
    void getInvPerm(int **invPerm, int *len);
    LevelData* getLevelDataRef();
    int getTotalLevel();
    int getTotalNodes();
    int* getLevelPtrRef();
    int* getNodePtrRef();
    int* getUnlockRowRef();
    int* getDangerRowRef();
    int* getUnlockCtrRef();

};

#endif
