/*
 * =======================================================================================
 *
 *   RACE: Recursicve Algebraic Coloring Engine
 *   Copyright (C) 2019, RRZE, Friedrich-Alexander-Universität Erlangen-Nürnberg
 *   Author: Christie Alappat
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * =======================================================================================
 */

#include "zone_tree.h"
#include "levelData.h"
#include "lb.h"
#include "utility.h"
#include <limits>

ZoneLeaf::ZoneLeaf():valueZ(2),nthreadsZ(-1),idealNthreadsZ(-1),totalSubBlocks(-1), parentZ(-2),effRowZ(-1),reachedLimit(false), stage(0), time(0)
{
}

ZoneLeaf::ZoneLeaf(int rangeLo_, int rangeHi_, int parent_):valueZ(2),nthreadsZ(1),idealNthreadsZ(1),totalSubBlocks(1), parentZ(parent_),  effRowZ(rangeHi_-rangeLo_),pinOrder(-1),reachedLimit(false), time(0)
{
    valueZ[0] = rangeLo_;
    valueZ[1] = rangeHi_;
}

ZoneTree::ZoneTree(RACE::dist dist_, RACE::d2Method d2Type_, RACE::LBTarget lbTarget_):cachedTree(NULL),dist(dist_),d2Type(d2Type_),lbTarget(lbTarget_),tree(NULL), maxStages_store(1)
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
    std::vector<int>* children = &(tree->at(parentIdx).children);
    //only till where children have been spawned
    int totalSubBlock = children->size()/2;
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
            for(int i=(children->at(2*subBlock)+block); i<children->at(2*subBlock+1); i+=blockPerThread)
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
    std::vector<int>* children = &(tree->at(parentIdx).children);

    //Only till where children have been spawned
    int totalSubBlocks = (children->size()/2);

    int nThreads = 0;
    int blockPerThread = getBlockPerThread(dist, d2Type);

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        //threads of subBlock
        int subNThreads = 0;
        int numThreads = (children->at(2*subBlock+1)-children->at(2*subBlock))/blockPerThread;//TODO for 3 block case

        for(int i=0; i<numThreads; ++i)
        {
            int maxNThreads = 0;
            for(int block=0; block<blockPerThread; ++block)
            {
                maxNThreads = std::max(cachedAt(children->at(2*subBlock)+blockPerThread*i+block).nthreadsZ, maxNThreads);
            }
            subNThreads += maxNThreads;
        }
        nThreads = std::max(nThreads, subNThreads);
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
    std::vector<int>* children = &(tree->at(currChild).children);

    int blockPerThread = getBlockPerThread(dist, d2Type);
    int totalSubBlocks = tree->at(currChild).totalSubBlocks;
    double *maxTime = (double*) malloc(sizeof(double)*blockPerThread);
    double sumTime = 0;

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
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
                int nxtChild = children->at(2*subBlock)+blockPerThread*i+block;
                updateTimeRecursive(nxtChild);
                maxTime[block] = std::max(maxTime[block], cachedAt(children->at(2*subBlock)+blockPerThread*i+block).time);
            }
        }
        sumTime += sumArr(maxTime, blockPerThread);
    }

    tree->at(currChild).time += sumTime;
    free(maxTime);
}

bool ZoneTree::spawnChild(int parentIdx, int parentSubIdx, int requestNthreads, LevelData* levelData, double eff)
{
//    int bestNThreads = 1;
//    int minEffRow = std::numeric_limits<int>::max();
    //TODO for three block
//    int maxThreads = 1;
    LB lb(requestNthreads, eff, levelData, dist, d2Type, lbTarget);
    lb.balance();

    int baseLen = tree->size();
    int *zonePtr = NULL, *subZonePtr = NULL, *numBlocks = NULL;
    int totalBlocks, totalSubBlocks;
    int base = tree->at(parentIdx).valueZ[parentSubIdx];
    int parent_stage = tree->at(parentIdx).stage;


    int len;
    lb.getZonePtr(&zonePtr, &len, base);
    //This gets a blocked version of zonePtr
    lb.getSubZonePtr(&subZonePtr, &totalSubBlocks, base);
    totalSubBlocks -= 1; //Since zonePtr would have one extra
    lb.getNumBlocks(&numBlocks, &totalBlocks);

    int ctr=0;

    tree->at(parentIdx).children.push_back(baseLen);

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
        currLeaf.stage = parent_stage+1;

        tree->push_back(currLeaf);
    }

    tree->at(parentIdx).children.push_back(baseLen+totalBlocks);


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
        printf("%d. Range:[", (int)i);

        for(unsigned j=0; j<currLeaf->valueZ.size(); ++j)
        {
            printf("%d, ",currLeaf->valueZ[j]);
        }

        printf("] Children:[");

        for(unsigned j=0; j<currLeaf->children.size(); ++j)
        {
            printf("%d, ",currLeaf->children[j]);
        }

        printf("], Parent: %d, nthreads: %d, idealNThreads:%d, stage:%d, subBLocks:%d, effRow: %d, reachedLimit: %s, pinOrder = %d funTime = %f\n", currLeaf->parentZ, currLeaf->nthreadsZ, currLeaf->idealNthreadsZ, currLeaf->stage, currLeaf->totalSubBlocks, currLeaf->effRowZ, (currLeaf->reachedLimit)?"true":"false", currLeaf->pinOrder, currLeaf->time);

        /*
        printf(", Pinned Core:[");
        for(unsigned j=0; j<currLeaf->pinnedCore.size(); ++j)
        {
            printf("%d, ",currLeaf->pinnedCore[j]);
        }
        printf("]\n");
        */

    }
}

void ZoneTree::findMaxStage()
{
    int maxStage = 0;
    for(unsigned i=0; i<tree->size(); ++i)
    {
        ZoneLeaf* currLeaf = &(at(i));
        maxStage = std::max(maxStage, currLeaf->stage);
    }

    //store it
    maxStages_store = maxStage;
}

int ZoneTree::maxStages()
{
    return maxStages_store;
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
