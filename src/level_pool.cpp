#include "level_pool.h"
#include "functional"

LevelPool::LevelPool(ZoneTree *zoneTree_, int SMT_, RACE::PinMethod pinMethod_):zoneTree(zoneTree_),pin(zoneTree_,SMT_, pinMethod_),tree(NULL)
{
    int totThreads = zoneTree->at(0).nthreadsZ;
    pool.init(totThreads);

    int ctr = 0;
    //create mappedIdx
    for(int i=0; i<zoneTree->size(); ++i)
    {
        mappedIdx.push_back(ctr);
        ctr += zoneTree->at(i).totalSubBlocks;
    }

    tree = new team<int> [ctr];
    pin.pinInit();
}

LevelPool::~LevelPool()
{
    delete[] tree;
}

//Recursively spawn thread pools
void LevelPool::createPoolRecursive(int parentIdx)
{
    RECURSIVE_HELPER(LevelPool::createPoolRecursive,
            std::vector<int> gid(nThreads);
            for(unsigned i=0; i<gid.size(); ++i)
            {
                gid[i] = zoneTree->at(children->at(2*parentSubIdx)+blockPerThread*i).pinOrder;
            }
            tree[poolTreeIdx(parentIdx, parentSubIdx)].init(gid, &pool);
    );
}

void LevelPool::createPool()
{
    int root = 0;
    createPoolRecursive(root);
    pinPool();
}


//Recursively pin pools
void LevelPool::pinPoolRecursive(int parentIdx)
{
    RECURSIVE_HELPER(LevelPool::pinPoolRecursive, pin.pinThread(zoneTree->at(parentIdx).pinOrder));
}


void LevelPool::pinPool()
{
    int root = 0;
    pinPoolRecursive(root);
}

void LevelPool::resetMaster()
{
    pin.resetMaster();
}

void LevelPool::sleepPoolRecursive(int parentIdx)
{
    RECURSIVE_HELPER(LevelPool::sleepPoolRecursive,
           tree[poolTreeIdx(parentIdx,parentSubIdx)].sleep();
    );
}


void LevelPool::sleepPool()
{
    int root = 0;
    sleepPoolRecursive(root);
}


