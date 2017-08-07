#include "pin.h"
#include "zone_tree.h"
#include "macros.h"

Pin::Pin(ZoneTree* zoneTree_, int SMT_, PinMethod method_):zoneTree(zoneTree_),machine(SMT_),SMT(SMT_),method(method_)
{

}

/*
void Pin::pinOrderRecursive(int parentIdx)
{
    ZoneLeaf* parentNode = &(zoneTree->at(parentIdx));
    std::vector<int> *children = &(parentNode->childrenZ);

    //TODO for 3 block scheme
    for(int block = 0; block<2; ++block)
    {
        int pinCtr = parentNode->pinOrder;
        int siblingInc = (block==0)?1:-1;
        for(unsigned i=block; i<children->size(); i+=2)
        {
            ZoneLeaf* childNode = &(zoneTree->at(children->at(i)));
            ZoneLeaf* siblingNode = &(zoneTree->at(children->at(i+siblingInc)));
            childNode->pinOrder = pinCtr;
            //Taking max since it might happen red and black have uneven threads
            pinCtr = pinCtr + std::max(childNode->nthreadsZ,siblingNode->nthreadsZ);
            pinOrderRecursive(children->at(i));
        }
    }
}
*/

void Pin::pinOrderRecursive(int parentIdx)
{
    ZoneLeaf* parentNode = &(zoneTree->at(parentIdx));
    std::vector<int> *children = &(parentNode->childrenZ);

    int blockPerThread = getBlockPerThread(zoneTree->dist, zoneTree->d2Type);
    int nThread = children->size()/blockPerThread;
    int pinCtr = parentNode->pinOrder;

    //TODO for 3 block scheme
    for(int i = 0; i<nThread; ++i)
    {
        int maxNThread = 0;

        //find max. of all siblings
        for(int block=0; block<blockPerThread; ++block)
        {
            ZoneLeaf* childNode = &(zoneTree->at(children->at(blockPerThread*i+block)));
            maxNThread = std::max(maxNThread, childNode->nthreadsZ);
        }

        for(int block=0; block<blockPerThread; ++block)
        {
            ZoneLeaf* childNode = &(zoneTree->at(children->at(blockPerThread*i+block)));
            childNode->pinOrder = pinCtr;
            pinOrderRecursive(children->at(blockPerThread*i+block));
        }
        pinCtr =  pinCtr + maxNThread;
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
    std::vector<int> *children = &(zoneTree->at(parent).childrenZ);

    int blockPerThread = getBlockPerThread(zoneTree->dist, zoneTree->d2Type);
    int currNthreads = static_cast<int>(children->size()/blockPerThread);

    if(currNthreads <= 1)
    {
        int tid = omp_get_thread_num();
        int pinOrder = zoneTree->at(parent).pinOrder;
#pragma omp critical
        {
        printf("%d ,Pinning: tid=%d -> cpu=%d",parent,tid,pinOrder);
        pinThread(pinOrder);
        printf("pinned to %d\n",sched_getcpu());
        }
    }
    else
    {
#pragma omp parallel num_threads(currNthreads)
        {
            int tid = omp_get_thread_num();
/*#pragma omp critical
            {
                printf("Pinning: tid=%d -> cpu=%d",tid,pinOrder);
                pinThread(pinOrder);
                printf("pinned to %d\n",sched_getcpu());
            }*/
            for(int block=0; block<blockPerThread; ++block)
            {
                pinApplicationRecursive(children->at(2*tid+block));
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
