#include "pin.h"
#include "zone_tree.h"
#include "macros.h"

Pin::Pin(ZoneTree* zoneTree_, int SMT_, PinMethod method_):zoneTree(zoneTree_),machine(SMT_),SMT(SMT_),method(method_)
{

}

void Pin::pinOrderRecursive(int parentIdx)
{
    ZoneLeaf* parentNode = &(zoneTree->at(parentIdx));
    std::vector<int> *children = &(parentNode->children);

    int blockPerThread = getBlockPerThread(zoneTree->dist, zoneTree->d2Type);
    int totalSubBlocks = parentNode->totalSubBlocks;

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        int pinCtr = parentNode->pinOrder;

        int nThread;
        if(!children->empty())
        {
            nThread = (children->at(2*subBlock+1) - children->at(2*subBlock))/blockPerThread;
        }
        else
        {
            nThread = 0;
        }

        for(int i = 0; i<nThread; ++i)
        {
            int maxNThread = 0;
            //find max. of all siblings
            for(int block=0; block<blockPerThread; ++block)
            {
                ZoneLeaf* childNode = &(zoneTree->at(children->at(2*subBlock)+blockPerThread*i+block));
                maxNThread = std::max(maxNThread, childNode->nthreadsZ);
            }

            for(int block=0; block<blockPerThread; ++block)
            {
                ZoneLeaf* childNode = &(zoneTree->at(children->at(2*subBlock)+blockPerThread*i+block));
                childNode->pinOrder = pinCtr;
                pinOrderRecursive(children->at(2*subBlock) + blockPerThread*i + block);
            }

            pinCtr =  pinCtr + maxNThread;
        }
    }
}


void Pin::calcPinOrder()
{
    int root = 0;
    (zoneTree->at(root)).pinOrder = 0;
    pinOrderRecursive(root);
}

void Pin::createPuNodeMapping()
{
    int totalThreads = zoneTree->at(0).nthreadsZ;
    int numNode = machine.getNumNode();

    std::vector<int> puPerNode(numNode,0);

    if(method == SCATTER)
    {
        //Assuming equal weight, TODO unequal weight
        int currPuPerNode = static_cast<int>(static_cast<double>(totalThreads)/numNode);

        int restPU = totalThreads%numNode;
        bool insufficientResource = false;
        for(int i=0;  i<restPU; ++i)
        {
            puPerNode[i] = currPuPerNode+1;
            if(puPerNode[i] > machine.getNumPuInNode(i))
            {
                insufficientResource = true;
            }

        }

        for(int i=restPU;  i<numNode; ++i)
        {
            puPerNode[i] = currPuPerNode;
            if(puPerNode[i] > machine.getNumPuInNode(i))
            {
                insufficientResource = true;
            }
        }

        if(insufficientResource)
        {
            ERROR_PRINT("You do not have sufficient resources to satisfy the threads");
        }
    }
    else
    {
        int restThread = totalThreads;
        int nodeId = 0;

        while( (restThread != 0) && (nodeId != numNode) )
        {
            int currPuPerNode = std::min(machine.getNumPuInNode(nodeId), restThread);

            puPerNode[nodeId] = currPuPerNode;
            restThread = (restThread - currPuPerNode);
            nodeId++;
        }

        if(restThread != 0)
        {
            ERROR_PRINT("You do not have sufficient resources to satisfy the threads");
        }
    }

    //Now create the map
    for(int i=0; i<numNode; ++i)
    {
        for(int j=0; j<puPerNode[i]; ++j)
        {
            std::pair<int,int> map(j,i);
            pinMap.push_back(map);
        }
    }

    if(totalThreads != static_cast<int>(pinMap.size()))
    {
        ERROR_PRINT("ERROR occured while creating thread pin map");
    }
}

void Pin::pinInit()
{
    calcPinOrder();
    createPuNodeMapping();
}

void Pin::pinThread(int pinOrder)
{
    int puId = pinMap[pinOrder].first;
    int nodeId = pinMap[pinOrder].second;

    machine.pinThread(puId,nodeId);
}

void Pin::resetMaster()
{
    machine.resetMaster();
}

//OMP nested pinning will not work
//just by pinning for the first time
//This function is deprecated
//refer pthreads pinning in level_pool.cpp
void Pin::pinApplicationRecursive(int parent)
{
    std::vector<int> *children = &(zoneTree->at(parent).children);
    int totalSubBlocks = zoneTree->at(parent).totalSubBlocks;
    int blockPerThread = getBlockPerThread(zoneTree->dist, zoneTree->d2Type);

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        int currNthreads = (children->at(2*subBlock+1) - children->at(2*subBlock))/blockPerThread;

        int tid = omp_get_thread_num();
        int pinOrder = zoneTree->at(parent).pinOrder;

        printf("%d ,Pinning: tid=%d -> cpu=%d",parent,tid,pinOrder);
        pinThread(pinOrder);
        printf("pinned to %d\n",sched_getcpu());


#pragma omp parallel num_threads(currNthreads)
        {
            int child_tid = omp_get_thread_num();
            for(int block=0; block<blockPerThread; ++block)
            {
                pinApplicationRecursive(children->at(2*subBlock)+2*child_tid+block);
#pragma omp barrier
            }
        }
    }
}

void Pin::pinApplication()
{
    calcPinOrder();
    int resetNestedState = omp_get_nested();
    int resetDynamicState = omp_get_dynamic();
    //set nested parallelism
    omp_set_nested(1);
    omp_set_dynamic(0);

    pinApplicationRecursive(0);

    omp_set_nested(resetNestedState);
    omp_set_dynamic(resetDynamicState);
}
