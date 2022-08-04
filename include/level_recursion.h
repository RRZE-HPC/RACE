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

#ifndef RACE_LEVEL_RECURSION_H
#define RACE_LEVEL_RECURSION_H
#include "graph.h"
#include "traverse.h"
#include "zone_tree.h"
#include "type.h"

class LevelRecursion
{
    private:
        RACE::Graph* graph;
        RACE::dist dist;
        RACE::d2Method d2Type;
        RACE::LBTarget lbTarget;
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
        LevelRecursion(RACE::Graph* graph_, int requestNThreads_, RACE::dist dist_, RACE::d2Method d2Type_=RACE::TWO_BLOCK, RACE::LBTarget lbTarget_=RACE::NNZ);
        ~LevelRecursion();
        void levelBalancing();
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getAvailableThreads();
        ZoneTree* getZoneTree();
};
#endif
