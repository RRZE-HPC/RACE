#include "level_pool.h"
#include "functional"

LevelPool::LevelPool(ZoneTree *zoneTree_, int SMT_, PinMethod pinMethod_):zoneTree(zoneTree_),pin(zoneTree_,SMT_, pinMethod_), tree(NULL)
{
    tree = new thpool<int> [zoneTree->size()];
    pin.pinInit();
}

LevelPool::~LevelPool()
{
    delete[] tree;
}


//Recursively spawn thread pools
void LevelPool::createPoolRecursive(int parentIdx)
{
    RECURSIVE_HELPER(LevelPool::createPoolRecursive, tree[parentIdx].init(nthreads));
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
