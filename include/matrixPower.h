#ifndef RACE_MTX_POWER_H
#define RACE_MTX_POWER_H

#include "print.h"
#include "error.h"
#include <vector>
#include "graph.h"
#include "traverse.h"
#include "levelData.h"

//this is per stage
class mtxPower
{
    Graph* graph;
    std::vector<int> levelPtr;
    std::vector<int> nodePtr;
    std::vector<int> unlockRow;//for p2p sync
    std::vector<int> dangerRow;// for p2p sync
    std::vector<int> unlockCtr;
    int* cacheLevelGroup;
    std::vector<int> hopelessRegions; //stores level index at hopeless boundary
    std::vector<int> hopelessNodePtr;
    std::vector<std::vector<int>> hopelessRegionPositiveBoundary; //stores rows in boundary
    std::vector<std::vector<int>> hopelessRegionNegativeBoundary; //stores rows in boundary
    LevelData* levelData;
    Traverse* traverser;
    int totalLevel;
    int highestPower;
    int numSharedCache;
    double cacheSize;
    double safetyFactor;
    int cache_violation_cutoff;
    int startRow;
    int endRow;

    void identifyHopelessRegions(std::vector<int> cacheViolatedLevel);
    void getHopelessStartEnd(int count, int *start, int *end);
    void getHopelessStartEnd(int count, int *start, int *end, std::vector<int> _hopelessRegions_);
    double getElemUpperLimit(int level);
    int workingBoundaryLength();

    public:
    mtxPower(Graph* graph_, int highestPower_, int numSharedCache, double cacheSize_, double safetyFactor_, int cache_violation_cutoff_, int startRow_, int endRow_);
    ~mtxPower();
    void findPartition();
    void splitSharedCacheDomain();
    std::vector<int> findMacroLevelPtr(int* zones);
    void consolidatePartition();
    void getStatNUMA();
    void findUnlockCtr();
    void createLevelPtr();
    double getBytePerNNZ();
    void powerRun(int power, int *rowPtr, int *col, double *A, double *x);
    void getPerm(int **perm, int *len);
    void getInvPerm(int **invPerm, int *len);
    LevelData* getLevelDataRef();
    int getTotalLevel();
    int getTotalNodes();
    std::vector<int> getLevelPtr();
    std::vector<int> getNodePtr();
    std::vector<int> getUnlockRow();
    std::vector<int> getDangerRow();
    std::vector<int> getUnlockCtr();
    std::vector<int> getHopelessRegions();
    std::vector<int> getHopelessNodePtr();
};


#endif
