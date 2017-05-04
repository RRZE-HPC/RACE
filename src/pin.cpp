#include "pin.h"

Pin::Pin(ZoneTree* zoneTree_, bool SMT_, PinMethod method_):zoneTree(zoneTree_),machine(SMT_),SMT(SMT_),method(method_)
{

}

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

    if(totalThreads != pinMap.size())
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

//OMP nested pinning will not working
//just by pinning for the first time
void Pin::pinApplicationRecursive(int parent)
{
    std::vector<int> *children = &(zoneTree->at(parent).childrenZ);
    int currNthreads = static_cast<int>(children->size()/2.0);

    if(currNthreads > 1)
    {
#pragma omp parallel num_threads(currNthreads)
        {
            int tid = omp_get_thread_num();
            int pinOrder = zoneTree->at(children->at(2*tid)).pinOrder;
#pragma omp critical
            {
                printf("Pinning: tid=%d -> cpu=%d",tid,pinOrder);
                pinThread(pinOrder);
                printf("pinned to %d\n",sched_getcpu());
            }
            pinApplicationRecursive(children->at(2*tid));
#pragma omp barrier
            pinApplicationRecursive(children->at(2*tid+1));
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
