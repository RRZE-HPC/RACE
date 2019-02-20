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
    int* macroLevelPtr;
    LevelData* levelData;
    Traverse* traverser;
    int totalLevel;
    int highestPower;
    double cacheSize;
    double safetyFactor;
    public:
    mtxPower(Graph* graph_, int highestPower_, double cacheSize_, double safetyFactor_);
    ~mtxPower();
    double getElemUpperLimit(int level);
    void findPartition();
    void consolidatePartition();
    void createLevelPtr();
    double getBytePerNNZ();
    void powerRun(int power, int *rowPtr, int *col, double *A, double *x);
    void getPerm(int **perm, int *len);
    void getInvPerm(int **invPerm, int *len);
    int getTotalLevel();
    LevelData* getLevelDataRef();
    int* getLevelPtrRef();
};

#endif
