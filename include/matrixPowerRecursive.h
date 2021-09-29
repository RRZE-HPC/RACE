#ifndef _MATRIX_POWER_RECURSIVE_H
#define _MATRIX_POWER_RECURSIVE_H

#include "matrixPower.h"
#include "string.h"

struct MPLeaf
{
    std::vector<int> nodePtr;
    int nodeId;

    //inputs
    int nrows;
    Range range;
    //boundaryRange structure
    //  ->workingBoundaryRadius (tells me what was its radius when generated
    //  initially from pure main region, this is necessary to limit the max power on
    //  these regions)
    //      ->radius (current radius) made a map, so can be empty too
    //          -> ranges
    //radii are stored -p -(p-1) ... -1 0 1 ... (p-1) (p)
    std::vector<std::map<int, std::vector<Range>>> boundaryRange;

    //outputs
    std::vector<int> hid; //hopelessRegions Ids;
    std::vector<int> lp;//levelPtr; //only of main
    //same structure as boundaryRange
    std::vector<std::map<int, std::vector<std::vector<int>>>> blp;
    std::vector<std::map<int, std::vector<std::vector<int>>>> boundaryUnlockRow;
    std::vector<std::map<int, std::vector<std::vector<int>>>> boundaryDangerRow;

    //std::vector<int> hopelessRegions;//I think this should be same as hid
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
    std::vector<int> childrenNodeStart;
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

    std::string mtxType;

    std::vector<int> cache_violation_cutoff;
    int get_cache_violation_cutoff(int stage);


    public:
    mtxPowerRecursive(Graph* graph_, int highestPower_, int numSharedCache, double cacheSize_, double safetyFactor_, std::string mtxType_="N");
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