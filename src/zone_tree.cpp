#include "zone_tree.h"
#include "levelData.h"
#include "lb.h"
#include <limits>

ZoneLeaf::ZoneLeaf():valueZ(2),nthreadsZ(-1),parentZ(-2),effRowZ(-1),reachedLimit(false)
{
}

ZoneLeaf::ZoneLeaf(int rangeLo_, int rangeHi_, int parent_):valueZ(2),nthreadsZ(1),parentZ(parent_),effRowZ(rangeHi_-rangeLo_),pinOrder(-1),reachedLimit(false)
{
    valueZ[0] = rangeLo_;
    valueZ[1] = rangeHi_;
}

ZoneTree::ZoneTree():cachedTree(NULL),tree(NULL)
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

KeyChild ZoneTree::findKeyChild(int parentIdx)
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
    keyChild.indices[0] = (children->at(2*maxIdx));
    keyChild.indices[1] = (children->at(2*maxIdx+1));
    keyChild.effRow = maxEffRow;
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
    int numThreads = static_cast<int>(children->size()/2.0);//TODO for 3 block case

    for(int i=0; i<numThreads; ++i)
    {
        int nThreadsRed = cachedAt(children->at(2*i)).nthreadsZ;
        int nThreadsBlack = cachedAt(children->at(2*i+1)).nthreadsZ;
        nThreads += std::max(nThreadsRed, nThreadsBlack);
    }
    tree->at(parentIdx).nthreadsZ = nThreads;

    int nxtParent = tree->at(parentIdx).parentZ;
    if(nxtParent!=-1)
    {
        updateTreeNThreads(nxtParent);
    }
}

bool ZoneTree::spawnChild(int parentIdx, int requestNthreads, int startThread, LevelData* levelData)
{
    int bestNThreads = 1;
    int minEffRow = std::numeric_limits<int>::max();
    //TODO for D1 and three block
    int maxThreads = static_cast<int>(levelData->totalLevel/4.0);
    int scanTill = std::min(requestNthreads, maxThreads);

    if( (startThread < maxThreads) && (scanTill > 1) )
    {
        tree->at(parentIdx).childrenZ.clear();
        ZoneLeaf resetParent = tree->at(parentIdx);

        ZoneLeaf bestParent;
        tree_t* bestParentSubTree = new tree_t;
        int baseLen = tree->size();

        for (int nthreads=startThread; nthreads<=scanTill; ++nthreads)
        {
            LB lb(nthreads, levelData);
            lb.D2LB();
            int *zonePtr = NULL;
            int len;
            int base = tree->at(parentIdx).valueZ[0];
            lb.getZonePtr(&zonePtr, &len, base);

            for(int block=0; block<2*nthreads; ++block)
            {
                ZoneLeaf currLeaf(zonePtr[block], zonePtr[block+1],parentIdx);

                cachedTree->push_back(currLeaf);
                tree->at(parentIdx).childrenZ.push_back(baseLen+block);
            }

            updateTreeEffRow(parentIdx);
            updateTreeNThreads(parentIdx);

            if(tree->at(parentIdx).effRowZ < minEffRow)
            {
                minEffRow = tree->at(0).effRowZ;
                bestNThreads = nthreads;
                //swap trees
                tree_t* tempTree;
                tempTree = bestParentSubTree;
                bestParentSubTree = cachedTree;
                cachedTree = tempTree;

                //correct parent
                bestParent = tree->at(parentIdx);
            }
            //clear parent
            tree->at(parentIdx) = resetParent;
            cachedTree->clear();
        }
        //insert bestSubTree into the main tree
        tree->at(parentIdx) = bestParent;
        tree->insert(tree->end(), bestParentSubTree->begin(), bestParentSubTree->end());
    }

    if(requestNthreads != bestNThreads)
    {
        return false;
    }
    else
    {
        return true;
    }
}

