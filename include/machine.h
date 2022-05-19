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

#ifndef RACE_MACHINE_H
#define RACE_MACHINE_H

#include <hwloc.h>
#include "error.h"
#include <vector>

struct PUNode{
    hwloc_cpuset_t puSet;
    int rank;
};

inline bool operator<(const PUNode &a, const PUNode &b)
{
    return (a.rank < b.rank);
}

class Machine{
    private:
        hwloc_topology_t topology;
        std::vector<std::vector<PUNode>> topTree;
        //required to reset the affinity
        hwloc_cpuset_t master_affinity;

        std::vector<int> corePerNode;
        int SMT;
        int numNode;
        int numCore;
        int numPU;
        void initTopology();
        void sortSMT();
    public:
        //constructor
        Machine(int SMT_);
        //Machine();
        ~Machine();
        int getNumNode();
        int getNumPU();
        int getNumCore();
        int getNumPuInNode(int logicalNodeid);
        int getNumCoreInNode(int logicalNodeid);
        void resetMaster();
        RACE_error pinThread(int logicalPUid, int logicalNodeId);
        void printTopTree();
};

#endif
