#include "pin.h"
#include "zone_tree.h"
#include "macros.h"
#include "machine.h"

Pin::Pin(ZoneTree* zoneTree_, int SMT_, RACE::PinMethod method_):zoneTree(zoneTree_),machine(NULL),SMT(SMT_),method(method_)
{
    Machine *mc = new Machine(SMT);
    machine = (void*) mc;
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

RACE_error Pin::createPuNodeMapping()
{
    Machine* mc = (Machine*) machine;
    int totalThreads = zoneTree->at(0).nthreadsZ;
    int numNode = mc->getNumNode();

    std::vector<int> puPerNode(numNode,0);

    if(method == RACE::SCATTER)
    {
        //Assuming equal weight, TODO unequal weight
        int currPuPerNode = static_cast<int>(static_cast<double>(totalThreads)/numNode);

        int restPU = totalThreads%numNode;
        bool insufficientResource = false;
        for(int i=0;  i<restPU; ++i)
        {
            puPerNode[i] = currPuPerNode+1;
            if(puPerNode[i] > mc->getNumPuInNode(i))
            {
                insufficientResource = true;
            }

        }

        for(int i=restPU;  i<numNode; ++i)
        {
            puPerNode[i] = currPuPerNode;
            if(puPerNode[i] > mc->getNumPuInNode(i))
            {
                insufficientResource = true;
            }
        }

        if(insufficientResource)
        {
            ERROR_PRINT("You do not have sufficient resources to satisfy the threads");
            return RACE_ERR_PIN;
        }
    }
    else
    {
        int restThread = totalThreads;
        int nodeId = 0;

        while( (restThread != 0) && (nodeId != numNode) )
        {
            int currPuPerNode = std::min(mc->getNumPuInNode(nodeId), restThread);

            puPerNode[nodeId] = currPuPerNode;
            restThread = (restThread - currPuPerNode);
            nodeId++;
        }

        if(restThread != 0)
        {
            ERROR_PRINT("You do not have sufficient resources to satisfy the threads");
            return RACE_ERR_PIN;
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
        return RACE_ERR_PIN;
    }

    return RACE_SUCCESS;
}

RACE_error Pin::pinInit()
{
    calcPinOrder();
    return createPuNodeMapping();
}

RACE_error Pin::pinThread(int pinOrder)
{
    Machine* mc = (Machine*) machine;
    int puId = pinMap[pinOrder].first;
    int nodeId = pinMap[pinOrder].second;

    return mc->pinThread(puId,nodeId);
}

void Pin::resetMaster()
{
    Machine* mc = (Machine*) machine;
    mc->resetMaster();
}

//OMP nested pinning will not work
//just by pinning for the first time
//This function is deprecated
//refer pthreads pinning in level_pool.cpp
RACE_error Pin::pinApplicationRecursive(int parent)
{
    std::vector<int> *children = &(zoneTree->at(parent).children);
    int totalSubBlocks = zoneTree->at(parent).totalSubBlocks;
    int blockPerThread = getBlockPerThread(zoneTree->dist, zoneTree->d2Type);
    RACE_error save_ret = RACE_SUCCESS;

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        int currNthreads;

        if(!children->empty())
        {
            currNthreads = (children->at(2*subBlock+1) - children->at(2*subBlock))/blockPerThread;
        }
        else
        {
            currNthreads = 0;
        }

        int pinOrder = zoneTree->at(parent).pinOrder;

        RACE_error ret = pinThread(pinOrder);
        if(ret != RACE_SUCCESS)
        {
            save_ret = ret;
        }

        if(currNthreads > 1)
        {
#pragma omp parallel num_threads(currNthreads)
            {
                int child_tid = omp_get_thread_num();
                for(int block=0; block<blockPerThread; ++block)
                {
                    ret = pinApplicationRecursive(children->at(2*subBlock)+blockPerThread*child_tid+block);
                    if(ret != RACE_SUCCESS)
                    {
                        save_ret = ret;
                    }
#pragma omp barrier
                }
            }
        }
    }

    return save_ret;
}

RACE_error Pin::pinApplication()
{
    calcPinOrder();
    int resetNestedState = omp_get_max_active_levels();
    int resetDynamicState = omp_get_dynamic();
    //set nested parallelism
    omp_set_max_active_levels(zoneTree->maxStages()+2); //+2 for safety
    omp_set_dynamic(0);

    RACE_error ret = pinApplicationRecursive(0);

    omp_set_max_active_levels(resetNestedState);
    omp_set_dynamic(resetDynamicState);

    return ret;
}

void Pin::pinPowerThread(int nodes)
{
/*    int resetNestedState = omp_get_max_active_levels();
    int resetDynamicState = omp_get_dynamic();
    //set nested parallelism
    //printf("setting nested\n");
    omp_set_max_active_levels(zoneTree->maxStages()+2);
    omp_set_dynamic(0);
    Machine* mc = (Machine*) machine;
    //body
#pragma omp parallel num_threads(nodes)
    {
        int parentId = omp_get_thread_num();
        mc->pinThread(0, parentId);
#pragma omp parallel
        {
            mc->pinThread(omp_get_thread_num(), parentId);
        }
    }

    //reset states
    omp_set_max_active_levels(resetNestedState);
    omp_set_dynamic(resetDynamicState);
    */
    UNUSED(nodes);
}

