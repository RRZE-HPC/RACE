#ifndef NAME_LEVEL_RECURSION_H
#define NAME_LEVEL_RECURSION_H
#include "graph.h"
#include "traverse.h"
#include "zone_tree.h"
#include "type.h"

class LevelRecursion
{
    private:
        Graph* graph;
        dist_t dist;
        d2Method d2Type;
        int requestNThreads;
        int availableNThreads;
        ZoneTree* zoneTree;
        int* perm;
        int* invPerm;
        void calculateIdealNthreads(int parentIdx, int currLvl);
        void recursivePartition(int parentIdx, int currLevel);
        //stores a vector of efficiency at different levels
        std::vector<double> eff_vec;
        std::vector<int> lvl_threads;
        double efficiency(unsigned levelNum);
        int lvlThreads(unsigned levelNum);
    public:
        LevelRecursion(Graph* graph_, int requestNThreads_, dist_t dist_, d2Method d2Type_=TWO_BLOCK);
        ~LevelRecursion();
        void levelBalancing();
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getAvailableThreads();
        ZoneTree* getZoneTree();
};
#endif
