#include "machine.h"
#include <algorithm>

Machine::Machine(int SMT_):SMT(SMT_)
{
    initTopology();
}

Machine::~Machine()
{
    hwloc_topology_destroy(topology);
}

void Machine::initTopology()
{
    hwloc_topology_init(&topology);
    hwloc_topology_load(topology);
    numNode = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PACKAGE);
    numNode = std::max(1,numNode);
    topTree.resize(numNode);
    corePerNode.resize(numNode);

    numCore = 0;
    numPU = 0;

   for(int i=0; i<numNode; ++i)
    {
        hwloc_cpuset_t totalCoreSet;

        totalCoreSet = hwloc_get_obj_by_type(topology, HWLOC_OBJ_PACKAGE, i)->cpuset;
        //SMT considered
        int coreInNode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, totalCoreSet, HWLOC_OBJ_CORE);
        numCore += coreInNode;
        corePerNode[i] = coreInNode;

        for(int j=0; j<coreInNode; ++j)
        {
            hwloc_cpuset_t currCpuSet = (hwloc_get_obj_inside_cpuset_by_type(topology, totalCoreSet, HWLOC_OBJ_CORE, j)->cpuset);

            int puInNode = 0;

            puInNode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, currCpuSet, HWLOC_OBJ_PU);
            puInNode = std::min(SMT,puInNode);

            numPU += puInNode;

            for(int k=0; k<puInNode; ++k)
            {
                PUNode puNode;
                puNode.rank = k;
                puNode.puSet = hwloc_get_obj_inside_cpuset_by_type(topology, currCpuSet, HWLOC_OBJ_PU, k)->cpuset;
                topTree[i].push_back(puNode);
            }
        }
    }

    if(SMT)
    {
        sortSMT();
    }

//    hwloc_get_cpubind(topology, master_affinity, HWLOC_CPUBIND_THREAD);
//    INFO_PRINT("Node = %d Core = %d PU = %d", numNode, numCore, numPU);
}

//Sorts PU within a node such that 
//first all PUsets having rank 0 appear, then rank 1 ....
//If one limits to rank 0; SMT is disabled
void Machine::sortSMT()
{
    for(int i=0; i<numNode; ++i)
    {
        std::stable_sort(topTree[i].begin(), topTree[i].end());
    }
}

int Machine::getNumNode()
{
    return numNode;
}

int Machine::getNumCore()
{
    return numCore;
}

int Machine::getNumPU()
{
    return numPU;
}

int Machine::getNumPuInNode(int logicalNodeId)
{
    return topTree[logicalNodeId].size();
}

int Machine::getNumCoreInNode(int logicalNodeId)
{
    return corePerNode[logicalNodeId];
}

//master should call it
void Machine::resetMaster()
{
    int err = hwloc_set_cpubind(topology, master_affinity, HWLOC_CPUBIND_THREAD);
    if(err == -1)
    {
        ERROR_PRINT("Master: thread binding could not be performed");
    }
}

//logicalCPUid: first cores numbered from 1:numCoreInNode then SMT
NAME_error Machine::pinThread(int logicalPUid, int logicalNodeid)
{
    NAME_error ret = NAME_SUCCESS;

    if( (logicalNodeid >= numNode) && (logicalPUid >= getNumPuInNode(logicalNodeid)) )
    {
        ERROR_PRINT("Node/PU requested to bind does not exist");
        ret = NAME_ERR_HWLOC;
        return ret;
    }

    hwloc_cpuset_t bindCPU = hwloc_bitmap_dup(topTree[logicalNodeid][logicalPUid].puSet);
    hwloc_bitmap_singlify(bindCPU);
    int err = hwloc_set_cpubind(topology, bindCPU, HWLOC_CPUBIND_THREAD);
    hwloc_bitmap_free(bindCPU);

   if(err == -1)
    {
        ERROR_PRINT("Thread binding could not be performed");
        ret = NAME_ERR_HWLOC;
    }

    return ret;
}

void Machine::printTopTree()
{
    printf("Num Node = %d, Num PU = %d\n",numNode, numPU);
    for(int i=0; i<numNode; ++i)
    {
        printf("%d , rank = ",i);
        for(int j=0; j<getNumPuInNode(i); ++j)
        {
            printf("%d, ", topTree[i][j].rank);
        }
        printf("\n");
    }
}
