#ifndef RACE_LEVEL_RECURSION_H
#define RACE_LEVEL_RECURSION_H
#include "graph.h"
#include "traverse.h"
#include "zone_tree.h"
#include "type.h"

class LevelRecursion
{
    private:
        Graph* graph;
        RACE::dist dist;
        RACE::d2Method d2Type;
        int requestNThreads;
        int availableNThreads;
        ZoneTree* zoneTree;
        int* perm;
        int* invPerm;
        void calculateIdealNthreads(int parentIdx, int parentSubIdx, int currLvl);
        void recursivePartition(int parentIdx, int parentSubIdx, int currLevel);
        //stores a vector of efficiency at different levels
        std::vector<double> eff_vec;
        std::vector<int> lvl_threads;
        double efficiency(unsigned levelNum);
        int lvlThreads(unsigned levelNum);
    public:
        LevelRecursion(Graph* graph_, int requestNThreads_, RACE::dist dist_, RACE::d2Method d2Type_=RACE::TWO_BLOCK);
        ~LevelRecursion();
        void levelBalancing();
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getAvailableThreads();
        ZoneTree* getZoneTree();
};
#endif
