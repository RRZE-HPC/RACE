#include "zone_tree.h"
#include "levelData.h"
#include "lb.h"
#include <limits>

ZoneLeaf::ZoneLeaf():valueZ(2),nthreadsZ(-1),idealNthreadsZ(-1),parentZ(-2),effRowZ(-1),reachedLimit(false), time(0)
{
}

ZoneLeaf::ZoneLeaf(int rangeLo_, int rangeHi_, int parent_):valueZ(2),nthreadsZ(1),idealNthreadsZ(1),parentZ(parent_),effRowZ(rangeHi_-rangeLo_),pinOrder(-1),reachedLimit(false), time(0)
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

void ZoneTree::updateTimeRecursive(int currChild)
{
    std::vector<int>* children = &(tree->at(currChild).childrenZ);

    double maxRedTime = cachedAt(currChild).time/2.0;
    double maxBlackTime = cachedAt(currChild).time/2.0;
    int numThreads = static_cast<int>(children->size()/2.0);//TODO for 3 block case

    for(int i=0; i<numThreads; ++i)
    {
	int nxtRedChild = children->at(2*i);
	int nxtBlackChild = children->at(2*i+1);
	updateTimeRecursive(nxtRedChild);
 	updateTimeRecursive(nxtBlackChild);
        maxRedTime = std::max(maxRedTime, cachedAt(children->at(2*i)).time);
        maxBlackTime = std::max(maxBlackTime, cachedAt(children->at(2*i+1)).time);
    }
    tree->at(currChild).time = maxRedTime+maxBlackTime;
}

bool ZoneTree::spawnChild(int parentIdx, int requestNthreads, int startThread, LevelData* levelData, dist_t dist, LBMode mode, double eff)
{
    int bestNThreads = 1;
    int minEffRow = std::numeric_limits<int>::max();
    //TODO for D1 and three block
    int maxThreads = static_cast<int>(levelData->totalLevel/4.0);
    int scanTill = std::min(requestNthreads, maxThreads);
    printf("maxThreads = %d\n",maxThreads);
    if( (startThread < maxThreads) && (scanTill > 1) )
    {
        ZoneLeaf bestParent = tree->at(parentIdx);
        tree->at(parentIdx).childrenZ.clear();
        ZoneLeaf resetParent = tree->at(parentIdx);

        tree_t* bestParentSubTree = new tree_t;
        int baseLen = tree->size();

        for (int nthreads=startThread; nthreads<=scanTill; ++nthreads)
        {
            LB lb(nthreads, levelData, dist);
            lb.balance();
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

            printf("eff Row = %d\n", tree->at(parentIdx).effRowZ);

            if(tree->at(parentIdx).effRowZ < minEffRow)
            {
		std::vector<int>* range = &(tree->at(parentIdx).valueZ);
		int currNrow = range->at(1) - range->at(0);
		double effNThread = currNrow/static_cast<double>(tree->at(parentIdx).effRowZ);
		printf("effNThreads = %f\n",effNThread);
		if((mode != EFFICIENCY) || (effNThread >= (eff/100.0)*nthreads) )
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
	    }
            //clear parent
            tree->at(parentIdx) = resetParent;
            cachedTree->clear();
        }
	printf("chosen bestNThreads = %d\n",bestNThreads);
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
