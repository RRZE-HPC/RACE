#ifndef _MATRIX_POWER_RECURSIVE_H
#define _MATRIX_POWER_RECURSIVE_H

#include "matrixPower.h"

struct MPLeaf
{
    int nrows;
    std::vector<int> range;
    std::vector<int> nodePtr;
    std::vector<int> levelPtr;
    std::vector<int> hopelessRegions;
    int nodeId;
    std::vector<int> unlockCtr;
    std::vector<int> unlockRow;
    std::vector<int> dangerRow;
    std::vector<int> lockCtr;
    std::vector<int> lockTableCtr;
    std::vector<int> unitPtr; //unit contains levelPtrs till next hopelessRegion and a hopelessRegion
    std::vector<int> unitNodePtr; //unit contains levelPtrs till next hopelessRegion and a hopelessRegion
    int parent;
    int stage;
    std::vector<int> children;
};

class mtxPowerRecursive
{
    Graph* graph;

    //final values; all this via tree
/*    int* levelPtr;
    int* nodePtr;
    int* unlockRow;//for p2p sync
    int* dangerRow;// for p2p sync
    int* unlockCtr;
*/

    int totalRows;
//    int totalLevel;
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

    //give public access
    std::vector<MPLeaf> tree;
    std::vector<int> hopelessNodePtr;


    void recursivePartition(int parentIdx);
    void findPartition();
    //these are final values
    void getPerm(int **perm, int *len);
    void getInvPerm(int **invPerm, int *len);
    LevelData* getLevelDataRef();
    void printTree();
    //int getTotalLevel();
    //int getTotalNodes();
    /*int* getLevelPtrRef();
    int* getUnlockRowRef();
    int* getDangerRowRef();
    int* getUnlockCtrRef();
    */
    //std::vector<int> getNodePtr();
    //std::vector<int> getHopelessNodePtr();
};

#endif
