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
