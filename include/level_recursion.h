#ifndef NAME_LEVEL_RECURSION_H
#define NAME_LEVEL_RECURSION_H
#include "graph.h"
#include "traverse.h"
#include "zone_tree.h"

class LevelRecursion
{
    private:
        Graph* graph;
        int requestNThreads;
        int availableNThreads;
        ZoneTree* zoneTree;
        int* perm;
        int* invPerm;
        bool recursivePartition(int parentIdx, int subRequestNThreads);
    public:
        LevelRecursion(Graph* graph_, int requestNThreads_);
        ~LevelRecursion();
        void levelBalancing();
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getAvailableThreads();
        ZoneTree* getZoneTree();
};
#endif
