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
        int requestNThreads;
        int availableNThreads;
        ZoneTree* zoneTree;
        int* perm;
        int* invPerm;
        bool recursivePartition(int parentIdx, int subRequestNThreads, int currLevel);
	//stores a vector of efficiency at different levels
	std::vector<double> eff_vec;
	double efficiency(unsigned levelNum);
    public:
        LevelRecursion(Graph* graph_, int requestNThreads_, dist_t dist_);
        ~LevelRecursion();
        void levelBalancing();
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getAvailableThreads();
        ZoneTree* getZoneTree();
};
#endif
