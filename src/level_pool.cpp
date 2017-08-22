#include "level_pool.h"
#include "functional"

LevelPool::LevelPool(ZoneTree *zoneTree_, int SMT_, PinMethod pinMethod_):zoneTree(zoneTree_),pin(zoneTree_,SMT_, pinMethod_),tree(NULL)
{
    int totThreads = zoneTree->at(0).nthreadsZ;
    pool.init(totThreads);
    tree = new team<int> [zoneTree->size()];
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
                gid[i] = zoneTree->at(children->at(blockPerThread*i)).pinOrder;
            }
            tree[parentIdx].init(gid, &pool);
    );
}

void LevelPool::createPool()
{
    int root = 0;
    createPoolRecursive(root);
    printf("created pool\n");
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
           tree[parentIdx].sleep();
    );
}


void LevelPool::sleepPool()
{
    int root = 0;
    sleepPoolRecursive(root);
}


