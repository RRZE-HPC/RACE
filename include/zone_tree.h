#ifndef NAME_ZONE_TREE_H
#define NAME_ZONE_TREE_H

#include "levelData.h"
#include "type.h"

struct ZoneLeaf{
    std::vector<int> valueZ;
    std::vector<int> childrenZ;
    int nthreadsZ;
    int parentZ;
    int effRowZ;
    int pinOrder;
    bool reachedLimit;
    ZoneLeaf();
    ZoneLeaf(int rangeLo, int rangeHi, int parentIdx);
};

struct KeyChild{
    std::vector<int> indices;//Red, Black Now for time-being; TODO for 3 block
    int effRow;
    KeyChild():indices(2)
    {
    };
};


class ZoneTree{
    private:
        /**
         * @brief cachedTree enables efficient rewinding
         */
        tree_t *cachedTree;

        ZoneLeaf& cachedAt(unsigned idx);
        void updateTreeEffRow(int parentIdx);
        void updateTreeNThreads(int parentIdx);

    public:
        tree_t *tree;

        ZoneTree();
        ~ZoneTree();

        ZoneLeaf& at(unsigned idx);
        void push_back(ZoneLeaf& leaf);
        int size();
        /**
         * @brief spawns child by callling load balancing
         * @param[in] parentIdx The Leaf idx from which to spawn
         * @param[in] requestNthreads number of threads requested
         * (Its just a request it can happen the parentIdx can't
         * spawn the request, actual value stored in nthreadsZ of the Leaf)
         * @param[in] levelData ptr to LevelData that stores information of the
         * level to be partitioned
         * param[out] boolean indicating whether the request could be satisfied
         */
        bool spawnChild(int parentIdx, int requestNthreads, int startThread, LevelData* levelData, dist_t dist);
        KeyChild findKeyChild(int parentIdx);
};

#endif
