#ifndef RACE_MTX_POWER_H
#define RACE_MTX_POWER_H

#include "print.h"
#include "error.h"
#include <vector>
#include "graph.h"
#include "traverse.h"
#include "levelData.h"

class mtxPower
{
    Graph* graph;
    int* levelPtr;
    int* levelGroupPtr;
    int* unlockRow;//for p2p sync
    int* dangerRow;// for p2p sync
    int* unlockCtr;
    int *cacheLevelGroup;
    LevelData* levelData;
    Traverse* traverser;
    int totalLevel;
    int highestPower;
    int numSharedCache;
    double cacheSize;
    double safetyFactor;
    std::vector<int> cache_violation_cutoff;
    int get_cache_violation_cutoff(int stage);
    public:
    mtxPower(Graph* graph_, int highestPower_, int numSharedCache, double cacheSize_, double safetyFactor_);
    ~mtxPower();
    double getElemUpperLimit(int level);
    std::vector<int> findPartition();
    void splitSharedCacheDomain();
    void findMacroLevelPtr(int* zones, int* macroLevelPtr);
    void consolidatePartition(std::vector<int> hopelessRegion);
    void getStatNUMA();
    void findUnlockCtr();
    void createLevelPtr();
    double getBytePerNNZ();
    void powerRun(int power, int *rowPtr, int *col, double *A, double *x);
    void getPerm(int **perm, int *len);
    void getInvPerm(int **invPerm, int *len);
    LevelData* getLevelDataRef();
    int getTotalLevel();
    int getTotalLevelGroup();
    int* getLevelPtrRef();
    int* getLevelGroupPtrRef();
    int* getUnlockRowRef();
    int* getDangerRowRef();
    int* getUnlockCtrRef();
};

std::vector<int> identifyHopelessRegions(std::vector<int> cacheViolatedLevel);
void getHopelessStartEnd(std::vector<int> hopelessRegion, int count, int *start, int *end);

#endif
