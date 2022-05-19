/*
 * =======================================================================================
 *
 *   RACE: Recursicve Algebraic Coloring Engine
 *   Copyright (C) 2019, RRZE, Friedrich-Alexander-Universität Erlangen-Nürnberg
 *   Author: Christie Alappat
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * =======================================================================================
 */

#ifndef RACE_ZONE_TREE_H
#define RACE_ZONE_TREE_H

#include "levelData.h"
#include "type.h"
#include "macros.h"

struct ZoneLeaf{
    std::vector<int> valueZ;
    std::vector<int> children;
    int nthreadsZ;
    int idealNthreadsZ;
    //total number of cache-blocked subBlocks
    //the leaf has to handle
    int totalSubBlocks;
    int parentZ;
    int effRowZ;
    int pinOrder;
    bool reachedLimit;
    int stage;
    double time;
    ZoneLeaf();
    ZoneLeaf(int rangeLo, int rangeHi, int parentIdx);
};

struct KeyChild {
    std::vector<int> indices;//Red, Black Now for time-being; TODO for 3 block
    int effRow;
    KeyChild(int numIndex):indices(numIndex)
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
        void updateTimeRecursive(int parentIdx);
        int maxStages_store;//calculate once and store, don't recalculate everytime

    public:
        RACE::dist dist;
        RACE::d2Method d2Type;
        RACE::LBTarget lbTarget;
        tree_t *tree;

        ZoneTree(RACE::dist dist, RACE::d2Method d2Type, RACE::LBTarget lbTarget);
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
         * @param[in] eff The minimum efficiency with which the current Lvl has
         * to be created
         * param[out] boolean indicating whether the request could be satisfied
         */
        bool spawnChild(int parentIdx, int parentSubIdx, int requestNthreads, LevelData* levelData, double eff=0);
        KeyChild findKeyChild(int parentIdx);
        void printTree();
        //calculates maxStages, done internally only once
        void findMaxStage();
        //returns the calculated maxStages value
        int maxStages();
        void resetTime();
        void updateTime();
};

#endif
