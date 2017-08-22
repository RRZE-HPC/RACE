#include "zone_tree.h"
#include "levelData.h"
#include "lb.h"
#include "utility.h"
#include <limits>

ZoneLeaf::ZoneLeaf():valueZ(2),nthreadsZ(-1),idealNthreadsZ(-1),parentZ(-2),effRowZ(-1),reachedLimit(false), time(0)
{
}

ZoneLeaf::ZoneLeaf(int rangeLo_, int rangeHi_, int parent_):valueZ(2),nthreadsZ(1),idealNthreadsZ(1),parentZ(parent_),effRowZ(rangeHi_-rangeLo_),pinOrder(-1),reachedLimit(false), time(0)
{
    valueZ[0] = rangeLo_;
    valueZ[1] = rangeHi_;
}

ZoneTree::ZoneTree(dist_t dist_, d2Method d2Type_):cachedTree(NULL),dist(dist_),d2Type(d2Type_),tree(NULL)
{
    tree = new tree_t;
    cachedTree = new tree_t;
}

ZoneTree::~ZoneTree()
{
    if(tree) {
        delete tree;
    }

    if(cachedTree) {
        delete cachedTree;
    }
}

ZoneLeaf& ZoneTree::at(unsigned idx)
{
    return tree->at(idx);
}

void ZoneTree::push_back(ZoneLeaf& leaf)
{
    tree->push_back(leaf);
}

int ZoneTree::size()
{
    return tree->size();
}

ZoneLeaf& ZoneTree::cachedAt(unsigned idx)
{
    if(cachedTree->empty())
    {
        return tree->at(idx);
    }
    else
    {
        if(idx < tree->size())
        {
            return tree->at(idx);
        }
        else
        {
            return cachedTree->at(idx-tree->size());
        }
    }
}

//wrong eff. row
/*KeyChild ZoneTree::findKeyChild(int parentIdx)
  {
  int maxEffRow = 0;
  std::vector<int>* children = &(tree->at(parentIdx).childrenZ);
  int numThreads = static_cast<int>(children->size()/2.0);//TODO for 3 block case
  int maxIdx = 0;

  for(int i=0; i<numThreads; ++i)
  {
  int nrowRed = cachedAt(children->at(2*i)).effRowZ;
  int nrowBlack = cachedAt(children->at(2*i+1)).effRowZ;

  int nrowTot = nrowRed+nrowBlack;
  if(nrowTot > maxEffRow)
  {
  maxEffRow = nrowTot;
  maxIdx = i;
  }
  }

  KeyChild keyChild;
//This is wrong since barrier has to be taken into account
keyChild.indices[0] = (children->at(2*maxIdx));
keyChild.indices[1] = (children->at(2*maxIdx+1));
keyChild.effRow = maxEffRow;
return keyChild;
}*/

/*
KeyChild ZoneTree::findKeyChild(int parentIdx)
{
    int maxEffRow = 0;
    int maxRowRed=0, maxRowBlack=0;
    std::vector<int>* children = &(tree->at(parentIdx).childrenZ);
    int numThreads = static_cast<int>(children->size()/2.0);//TODO for 3 block case
    int maxIdx = 0;
    int maxIdxRed = 0, maxIdxBlack = 0;

    for(int i=0; i<numThreads; ++i)
    {
        int nrowRed = cachedAt(children->at(2*i)).effRowZ;

        if(nrowRed > maxRowRed)
        {
            maxRowRed = nrowRed;
            maxIdxRed = i;
        }
    }

    for(int i=0; i<numThreads; ++i)
    {
        int nrowBlack = cachedAt(children->at(2*i+1)).effRowZ;

        if(nrowBlack > maxRowBlack)
        {
            maxRowBlack = nrowBlack;
            maxIdxBlack = i;
        }
    }

    maxEffRow = maxRowRed + maxRowBlack;
    maxIdx = (maxRowRed>=maxRowBlack)?maxIdxRed:maxIdxBlack;
    KeyChild keyChild(2);
    keyChild.indices[0] = (children->at(2*maxIdx));
    keyChild.indices[1] = (children->at(2*maxIdx+1));
    keyChild.effRow = maxEffRow;
    return keyChild;
}
*/

KeyChild ZoneTree::findKeyChild(int parentIdx)
{
    int maxEffRow = 0;
    std::vector<int>* children = &(tree->at(parentIdx).childrenZ);
    int blockPerThread = getBlockPerThread(dist, d2Type);
    int numThreads = static_cast<int>(children->size()/blockPerThread);//TODO for 3 block case
    int maxTotIdx = 0;
    int *maxRow = new int[blockPerThread];
    int *maxIdx = new int[blockPerThread];

    for(int block=0; block<blockPerThread; ++block)
    {
        maxRow[block] = 0;
        maxIdx[block] = 0;
    }

    for(int block=0; block<blockPerThread; ++block)
    {
        for(int i=0; i<numThreads; ++i)
        {
            int nrow = cachedAt(children->at(blockPerThread*i+block)).effRowZ;

            if(nrow > maxRow[block])
            {
                maxRow[block] = nrow;
                maxIdx[block] = i;
            }
        }
    }

    maxEffRow = sumArr(maxRow, blockPerThread);
    maxTotIdx = maxArr(maxIdx, blockPerThread);
    KeyChild keyChild(blockPerThread);
    for(int block=0; block<blockPerThread; ++block)
    {
        keyChild.indices[block] = (children->at(blockPerThread*maxTotIdx)+block);
    }

    keyChild.effRow = maxEffRow;
    delete[] maxRow;
    delete[] maxIdx;
    return keyChild;
}


void ZoneTree::updateTreeEffRow(int parentIdx)
{
    KeyChild keyChild = findKeyChild(parentIdx);
    tree->at(parentIdx).effRowZ = keyChild.effRow;

    int nxtParent = tree->at(parentIdx).parentZ;
    if(nxtParent!=-1)
    {
        updateTreeEffRow(nxtParent);
    }
}

void ZoneTree::updateTreeNThreads(int parentIdx)
{
    std::vector<int>* children = &(tree->at(parentIdx).childrenZ);

    int nThreads = 0;
    int blockPerThread = getBlockPerThread(dist, d2Type);
    int numThreads = static_cast<int>(children->size()/blockPerThread);//TODO for 3 block case

    for(int i=0; i<numThreads; ++i)
    {
       int maxNThreads = 0;
       for(int block=0; block<blockPerThread; ++block)
       {
           maxNThreads = std::max(cachedAt(children->at(blockPerThread*i+block)).nthreadsZ, maxNThreads);
       }
       nThreads += maxNThreads;
    }
    tree->at(parentIdx).nthreadsZ = nThreads;

    int nxtParent = tree->at(parentIdx).parentZ;
    if(nxtParent!=-1)
    {
        updateTreeNThreads(nxtParent);
    }
}

void ZoneTree::updateTimeRecursive(int currChild)
{
    std::vector<int>* children = &(tree->at(currChild).childrenZ);

    int blockPerThread = getBlockPerThread(dist, d2Type);
    double *maxTime = (double*) malloc(sizeof(double)*blockPerThread);
    for(int block=0; block<blockPerThread; ++block)
    {
        //this is the time if load is perfectly balanced
        maxTime[block] = cachedAt(currChild).time/((double)blockPerThread);
    }

    int numThreads = static_cast<int>(children->size()/blockPerThread);//TODO for 3 block case

    for(int i=0; i<numThreads; ++i)
    {
        for(int block=0; block<blockPerThread; ++block)
        {
            int nxtChild = children->at(blockPerThread*i+block);
            updateTimeRecursive(nxtChild);
            maxTime[block] = std::max(maxTime[block], cachedAt(children->at(blockPerThread*i+block)).time);
        }
    }
    tree->at(currChild).time += sumArr(maxTime, blockPerThread);
    free(maxTime);
}

bool ZoneTree::spawnChild(int parentIdx, int requestNthreads, LevelData* levelData, double eff)
{
//    int bestNThreads = 1;
//    int minEffRow = std::numeric_limits<int>::max();
    //TODO for three block
//    int maxThreads = 1;

    printf("requesting %d threads\n", requestNthreads);
    LB lb(requestNthreads, eff, levelData, dist, d2Type);
    lb.balance();

    int baseLen = tree->size();
    int *zonePtr = NULL;
    int len;
    int base = tree->at(parentIdx).valueZ[0];
    lb.getZonePtr(&zonePtr, &len, base);

    for(int block=0; block<len-1; ++block)
    {
        ZoneLeaf currLeaf(zonePtr[block], zonePtr[block+1], parentIdx);

        tree->push_back(currLeaf);
        tree->at(parentIdx).childrenZ.push_back(baseLen+block);
    }

    updateTreeEffRow(parentIdx);
    updateTreeNThreads(parentIdx);

    printf("eff Row = %d\n", tree->at(parentIdx).effRowZ);

    int currNThreads = lb.getCurrLvlNThreads();

    printf("currNThreads = %d\n", currNThreads);

    if(currNThreads != requestNthreads)
    {
        return false;
    }
    else
    {
        return true;
    }

}

void ZoneTree::printTree()
{
    updateTime();
    for(unsigned i=0; i<tree->size(); ++i)
    {
        ZoneLeaf* currLeaf = &(at(i));
        printf("%d. Range:[", i);

        for(unsigned j=0; j<currLeaf->valueZ.size(); ++j)
        {
            printf("%d, ",currLeaf->valueZ[j]);
        }

        printf("] Children:[");


        for(unsigned j=0; j<currLeaf->childrenZ.size(); ++j)
        {
            printf("%d, ",currLeaf->childrenZ[j]);
        }
        printf("], Parent: %d, nthreads: %d, idealNThreads:%d, effRow: %d reachedLimit: %s, pinOrder = %d funTime = %f\n", currLeaf->parentZ, currLeaf->nthreadsZ, currLeaf->idealNthreadsZ, currLeaf->effRowZ, (currLeaf->reachedLimit)?"true":"false", currLeaf->pinOrder, currLeaf->time);
    }
}


void ZoneTree::resetTime()
{
    for(unsigned i=0; i<tree->size(); ++i)
    {
        tree->at(i).time = 0;
    }
}

void ZoneTree::updateTime()
{
    updateTimeRecursive(0);
} 
