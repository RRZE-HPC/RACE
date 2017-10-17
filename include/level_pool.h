#ifndef RACE_LEVEL_POOL_H
#define RACE_LEVEL_POOL_H

#include "thpool.h"
#include "pin.h"
#include "type.h"

/*
#define RECURSIVE_HELPER(this_fn, ...) \
    std::vector<int>* children = &(zoneTree->at(parentIdx).childrenZ);\
    int blockPerThread=2;\
    int nThreads = children->size()/blockPerThread;\
    __VA_ARGS__;\
    if(nThreads > 1)\
    {\
        for(int tid=0; tid<nThreads; ++tid)\
        {\
            tree[parentIdx].addJob(tid, std::bind(&this_fn,this, std::placeholders::_1), children->at(2*tid));\
        }\
        tree[parentIdx].barrier();\
        for(int tid=0; tid<nThreads; ++tid)\
        {\
            tree[parentIdx].addJob(tid, std::bind(&this_fn,this, std::placeholders::_1), children->at(2*tid+1));\
        }\
        tree[parentIdx].barrier();\
    }\
*/

#define RECURSIVE_HELPER(this_fn, ...) \
    std::vector<int>* children = &(zoneTree->at(parentIdx).children);\
    int blockPerThread = getBlockPerThread(zoneTree->dist, zoneTree->d2Type);\
    int totalSubBlocks = zoneTree->at(parentIdx).totalSubBlocks;\
    for(int parentSubIdx=0; parentSubIdx<totalSubBlocks; ++parentSubIdx)\
    {\
        int nThreads;\
        if(!children->empty())\
        {\
            nThreads = (children->at(2*parentSubIdx+1) - children->at(2*parentSubIdx))/blockPerThread;\
        }\
        else\
        {\
            nThreads = 0;\
        }\
        __VA_ARGS__;\
        if(nThreads > 1)\
        {\
            for(int block=0; block<blockPerThread; ++block)\
            {\
                for(int tid=0; tid<nThreads; ++tid)\
                {\
                    tree[poolTreeIdx(parentIdx, parentSubIdx)].addJob(tid, std::bind(&this_fn,this, std::placeholders::_1), children->at(2*parentSubIdx) + blockPerThread*tid + block);\
                }\
                tree[poolTreeIdx(parentIdx, parentSubIdx)].barrier();\
            }\
        }\
   }\



/*
#define RECURSIVE_HELPER(this_fn, fn) \
    std::vector<int>* children = &(zoneTree->at(parentIdx).childrenZ);\
    int nthreads = children->size()/2;\
    fn;\
    if(nthreads > 1)\
    {\
        std::vector<int> args(nthreads);\
        for(int tid=0; tid<nthreads; ++tid)\
        {\
            args[tid] = children->at(2*tid);\
        }\
        tree[parentIdx].addJob(std::bind(&this_fn,this, std::placeholders::_1), args);\
        tree[parentIdx].doJob();\
        for(int tid=0; tid<nthreads; ++tid)\
        {\
            args[tid] = children->at(2*tid+1);\
        }\
        tree[parentIdx].addJob(std::bind(&this_fn,this, std::placeholders::_1), args);\
        tree[parentIdx].doJob();\
    }\
*/

class LevelPool{
    private:
        ZoneTree* zoneTree;
        void createPoolRecursive(int parentIdx);
        void pinPoolRecursive(int parentIdx);
        void sleepPoolRecursive(int parentIdx);
    public:
        LevelPool(ZoneTree *zoneTree_, int SMT, PinMethod pinMethod);
        ~LevelPool();
        Pin pin;
        //thread pool
        thpool<int> pool;
        //Team tree
        team<int>* tree;
        //to map between (parentIdx, parentSubIdx) to poolTreeIdx
        //see macro poolTreeIdx
        std::vector<int> mappedIdx;
        inline int poolTreeIdx(int parentIdx, int parentSubIdx)
        {
            return (mappedIdx[parentIdx]+parentSubIdx);
        }
        //creates pinned pool
        void createPool();
        void pinPool();
        void sleepPool();
        void wake();
        void resetMaster();
};

struct Arg{
    LevelPool *pool;
    int parentIdx;
};


#endif
