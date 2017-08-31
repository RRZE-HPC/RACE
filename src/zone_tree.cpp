#include "zone_tree.h"
#include "levelData.h"
#include "lb.h"
#include "utility.h"
#include <limits>

ZoneLeaf::ZoneLeaf():valueZ(2),nthreadsZ(-1),idealNthreadsZ(-1),totalSubBlocks(-1), parentZ(-2),effRowZ(-1),reachedLimit(false), time(0)
{
}

ZoneLeaf::ZoneLeaf(int rangeLo_, int rangeHi_, int parent_):valueZ(2),nthreadsZ(1),idealNthreadsZ(1),totalSubBlocks(1), parentZ(parent_),  effRowZ(rangeHi_-rangeLo_),pinOrder(-1),reachedLimit(false), time(0)
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
    std::vector<int>* subPointer = &(tree->at(parentIdx).subPointer);
    //only till where children have been spawned
    int totalSubBlock = subPointer->size()/2;
    int blockPerThread = getBlockPerThread(dist, d2Type);
    int maxTotIdx = 0;
    int *maxRow = new int[blockPerThread];
    int *maxIdx = new int[blockPerThread];

    for(int block=0; block<blockPerThread; ++block)
    {
        maxRow[block] = 0;
        maxIdx[block] = 0;
    }

    for(int subBlock=0; subBlock<totalSubBlock; ++subBlock)
    {
        int *tempMaxRow = new int[blockPerThread];
        int *tempMaxIdx = new int[blockPerThread];

        for(int block=0; block<blockPerThread; ++block)
        {
            tempMaxRow[block] = 0;
            tempMaxIdx[block] = 0;
        }

        for(int block=0; block<blockPerThread; ++block)
        {
            for(int i=(subPointer->at(2*subBlock)+block); i<subPointer->at(2*subBlock+1); i+=blockPerThread)
            {
                int nrow = cachedAt(i).effRowZ;

                if(nrow > tempMaxRow[block])
                {
                    tempMaxRow[block] = nrow;
                    tempMaxIdx[block] = i;
                }
            }
        }

        for(int block=0; block<blockPerThread; ++block)
        {
            maxRow[block] += tempMaxRow[block];
            maxIdx[block] = tempMaxIdx[block]; //TODO
        }

        delete[] tempMaxRow;
        delete[] tempMaxIdx;
    }

    printf("finished\n");
    maxEffRow = sumArr(maxRow, blockPerThread);
    maxTotIdx = maxArr(maxIdx, blockPerThread);
    KeyChild keyChild(blockPerThread);
    for(int block=0; block<blockPerThread; ++block)
    {
        keyChild.indices[block] = maxTotIdx;
    }

    keyChild.effRow = maxEffRow;
    delete[] maxRow;
    delete[] maxIdx;
    printf("finished\n");
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
    std::vector<int>* subPointer = &(tree->at(parentIdx).subPointer);

    //Only till where children have been spawned
    int totalSubBlocks = (subPointer->size()/2);

    int nThreads = 0;
    int blockPerThread = getBlockPerThread(dist, d2Type);

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        //threads of subBlock
        int subNThreads = 0;
        int numThreads = (subPointer->at(2*subBlock+1)-subPointer->at(2*subBlock))/blockPerThread;//TODO for 3 block case

        for(int i=0; i<numThreads; ++i)
        {
            int maxNThreads = 0;
            for(int block=0; block<blockPerThread; ++block)
            {
                maxNThreads = std::max(cachedAt(subPointer->at(2*subBlock)+blockPerThread*i+block).nthreadsZ, maxNThreads);
            }
            subNThreads += maxNThreads;
        }

        printf("subBlock = %d, numThreads = %d subNThreads = %d\n", subBlock, numThreads, subNThreads);

        nThreads = std::max(nThreads, subNThreads);
    }

    printf("ParentIdx = %d, nThreads = %d\n", parentIdx, nThreads);
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

bool ZoneTree::spawnChild(int parentIdx, int parentSubIdx, int requestNthreads, LevelData* levelData, double eff)
{
//    int bestNThreads = 1;
//    int minEffRow = std::numeric_limits<int>::max();
    //TODO for three block
//    int maxThreads = 1;
    printf("requesting %d threads\n", requestNthreads);
    LB lb(requestNthreads, eff, levelData, dist, d2Type);
    lb.balance();

    int baseLen = tree->size();
    int *zonePtr = NULL, *subZonePtr = NULL, *numBlocks = NULL;
    int totalBlocks, totalSubBlocks;
    int base = tree->at(parentIdx).valueZ[parentSubIdx];

    int len;
    lb.getZonePtr(&zonePtr, &len, base);
    //This gets a blocked version of zonePtr
    lb.getSubZonePtr(&subZonePtr, &totalSubBlocks, base);
    totalSubBlocks -= 1; //Since zonePtr would have one extra
    lb.getNumBlocks(&numBlocks, &totalBlocks);

    printf("totalBlocks = %d, totalSubBlocks = %d\n", totalBlocks, totalSubBlocks);
    int ctr=0;

    tree->at(parentIdx).subPointer.push_back(baseLen);

    for(int block=0; block<totalBlocks; ++block)
    {
        ZoneLeaf currLeaf(zonePtr[block], zonePtr[block+1], parentIdx);
        currLeaf.valueZ.resize(numBlocks[block]+1);
        for(int subBlock=0; subBlock<numBlocks[block]; ++subBlock)
        {
            currLeaf.valueZ[subBlock] = (subZonePtr[ctr]);
            ++ctr;
        }
        currLeaf.valueZ[numBlocks[block]] = (subZonePtr[ctr]);
        currLeaf.totalSubBlocks = numBlocks[block];

        tree->push_back(currLeaf);
        tree->at(parentIdx).childrenZ.push_back(baseLen+block);
    }

    tree->at(parentIdx).subPointer.push_back(baseLen+totalBlocks);


    if(parentSubIdx == (tree->at(parentIdx).totalSubBlocks-1))
    {
        updateTreeEffRow(parentIdx);
        updateTreeNThreads(parentIdx);
    }

    int currNThreads = lb.getCurrLvlNThreads();

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

        printf("] SubPointer:[");


        for(unsigned j=0; j<currLeaf->subPointer.size(); ++j)
        {
            printf("%d, ",currLeaf->subPointer[j]);
        }

        printf("], Parent: %d, nthreads: %d, idealNThreads:%d, subBLocks:%d, effRow: %d, reachedLimit: %s, pinOrder = %d funTime = %f\n", currLeaf->parentZ, currLeaf->nthreadsZ, currLeaf->idealNthreadsZ, currLeaf->totalSubBlocks, currLeaf->effRowZ, (currLeaf->reachedLimit)?"true":"false", currLeaf->pinOrder, currLeaf->time);
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
