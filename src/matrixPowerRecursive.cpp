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

#include "matrixPowerRecursive.h"
#include "utility.h"
#include "config.h"

#define UPDATE_INVPERM(_invPerm_, _perm_)\
{\
    _Pragma("omp parallel for schedule(static)")\
    for(int i=0; i<totalRows; ++i)\
    {\
        _invPerm_[_perm_[i]] = i;\
    }\
}\

mtxPowerRecursive::mtxPowerRecursive(RACE::Graph* graph_, int highestPower_, int numSharedCache_, double cacheSize_, double safetyFactor_, std::string mtxType_):graph(graph_), highestPower(highestPower_), numSharedCache(numSharedCache_), cacheSize(cacheSize_), safetyFactor(safetyFactor_), perm(NULL), invPerm(NULL), mtxType(mtxType_)
{
    std::vector<int> default_cutoff;
    getEnv("RACE_CACHE_VIOLATION_CUTOFF_DEFAULT", default_cutoff);
    if(default_cutoff.empty())
    {
        default_cutoff.push_back(50);//make 50 factor so not much cut-off happens by default
    }
    cache_violation_cutoff.clear();
    cache_violation_cutoff.push_back(default_cutoff[0]); //default
    getEnv("RACE_CACHE_VIOLATION_CUTOFF", cache_violation_cutoff);

    maxRecStages = 0; //0 => no recursion
    std::vector<int> maxRecStagesVec;
    getEnv("RACE_MAX_RECURSION_STAGES", maxRecStagesVec);
    if(!maxRecStagesVec.empty())
    {
        maxRecStages = maxRecStagesVec[0];
    }
    else
    {
        //one specified by RACE_CACHE_VIOLATION_CUTOFF
        maxRecStages = std::max(0, (int)(cache_violation_cutoff.size()-1+1)); //+1 because cache violation is specified for current stage and depending on the value one more recursion stage can be done
    }

    totalRows = graph->NROW;
#ifndef RACE_PERMUTE_ON_FLY
    perm =  new int[totalRows];
    invPerm =  new int[totalRows];
    //copy initial permutations
#pragma omp parallel for schedule(static)
    for(int i=0; i<totalRows; ++i)
    {
        perm[i] = graph->totalPerm[i];
        invPerm[i] = graph->totalInvPerm[i];
    }
#endif
}

mtxPowerRecursive::~mtxPowerRecursive()
{
    if(perm)
    {
        delete[] perm;
    }
    if(invPerm)
    {
        delete[] invPerm;
    }
}

int mtxPowerRecursive::get_cache_violation_cutoff(int stage)
{
    int cache_cutoff = 50;

    if(stage >= maxRecStages)
    {
        cache_cutoff = std::numeric_limits<int>::max(); //disable further recursion
    }
    else
    {
        if(stage < (int)(cache_violation_cutoff.size()-1))
        {
            cache_cutoff = cache_violation_cutoff[stage+1]; //Stage numbering starts with 0, therefore stage+1
        }
        else
        {
            cache_cutoff = cache_violation_cutoff[0];
        }
    }

    return cache_cutoff;
}

#define READ_STAGE_DATA(_leaf_, _stage_mtxPower_, _readPermute_)\
    _leaf_.hid = _stage_mtxPower_.getHopelessRegions();\
    _leaf_.lp = _stage_mtxPower_.getLevelPtr();\
    _leaf_.blp = _stage_mtxPower_.getBoundaryLevelPtr();\
    _leaf_.unlockRow = _stage_mtxPower_.getUnlockRow();\
    _leaf_.boundaryUnlockRow = _stage_mtxPower_.getBoundaryUnlockRow();\
    _leaf_.unlockCtr = _stage_mtxPower_.getUnlockCtr();\
    _leaf_.dangerRow = _stage_mtxPower_.getDangerRow();\
    _leaf_.boundaryDangerRow = _stage_mtxPower_.getBoundaryDangerRow();\
    if(_readPermute_)\
    {\
        _stage_mtxPower_.getPerm(&perm_curStage, &_leaf_.nrows);\
    }\
    std::vector<int> _nodePtr_ = _stage_mtxPower_.getNodePtr();\
    _leaf_.nodePtr = _nodePtr_;\
    /*create unitPtr*/\
    int _n_nodes_ = _nodePtr_.size()-1;\
    std::vector<int> _hopelessNodePtr_ = _stage_mtxPower_.getHopelessNodePtr();\
    _leaf_.unitNodePtr.push_back(0);\
    _leaf_.unitPtr.push_back(_nodePtr_[0]);\
    for(int n=0; n<_n_nodes_; ++n)\
    {\
        /*int levelNodeStart = _nodePtr_[n];*/\
        int levelNodeEnd = _nodePtr_[n+1];\
        for(int hn=_hopelessNodePtr_[n]; hn<_hopelessNodePtr_[n+1]; ++hn)\
        {\
            printf("@@@@@@@@@ hopelessNodePtr[%d] = %d, hopelessNodePtr[%d] = %d\n", n, _hopelessNodePtr_[n], n+1, _hopelessNodePtr_[n+1]);\
            _leaf_.unitPtr.push_back(_leaf_.hid[hn]); \
            _leaf_.unitPtr.push_back(_leaf_.hid[hn]+1); \
        }\
        if(_leaf_.unitPtr.back() != levelNodeEnd) /*if there is no hopeless in last unit, make a dummy*/\
        {\
            _leaf_.unitPtr.push_back(levelNodeEnd);\
            _leaf_.unitPtr.push_back(levelNodeEnd);\
        }\
        _leaf_.unitNodePtr.push_back(_leaf_.unitPtr.size()-1);\
    }

#define COPY_TO_PLAIN_INT_PTR(_vec_, _plain_)\
    _plain_ = new int[_vec_.size()];\
    for(int i=0; i<(int)_vec_.size(); ++i)\
    {\
        _plain_[i] = _vec_[i];\
    }\

void mtxPowerRecursive::recursivePartition(int parentIdx)
{
    MPLeaf parentLeaf = tree[parentIdx];
    std::vector<int> parentHopelessRegions = parentLeaf.hid;
    int n_hopeless = (int)(parentHopelessRegions.size());
    int parentStage = parentLeaf.stage;

    std::vector<int> parentNodePtr = parentLeaf.nodePtr;
    int totalNodes = (int)(parentNodePtr.size()-1);
    if(parentStage==0)
    {
        tree[parentIdx].childrenNodeStart.resize(totalNodes,0);
    }
    else
    {
        tree[parentIdx].childrenNodeStart.resize(1,0);
    }

    for(int h=0; h<n_hopeless; ++h)
    {
        MPLeaf curLeaf;
        curLeaf.parent = parentIdx;
        curLeaf.stage = parentStage + 1;
        int hopelessStart = parentHopelessRegions[h];
        printf("@@@ hopless = %d,%d\n", hopelessStart, hopelessStart+1);
        int push_nodeId = 0;
        if(parentStage>0)
        {
            curLeaf.nodeId = parentLeaf.nodeId;

        }
        else
        {
            //search and find
            for(int n=0; n<totalNodes; ++n)
            {
                if( (parentNodePtr[n] <= hopelessStart) && (parentNodePtr[n+1] > hopelessStart) )
                {
                    curLeaf.nodeId = n;
                    break;
                }
            }
            push_nodeId = curLeaf.nodeId;
        }
        curLeaf.range.lo = parentLeaf.lp[hopelessStart];
        curLeaf.range.hi = parentLeaf.lp[hopelessStart+1];
        int wbl = workingBoundaryLength_base(highestPower);

        curLeaf.boundaryRange.resize(wbl);
        if(!parentLeaf.blp.empty())
        {
            for(int wr=0; wr<wbl; ++wr)//working radius will store the max radius parent has attained, so to limit the powers computed
            {
                //std::map<int, std::vector<std::vector<int>>> *curParentBLP = curParentBLP = &(parentLeaf.blp[wr]);
                std::map<int, std::vector<std::vector<int>>> *curParentBLP = &(parentLeaf.blp[wr]);
                for(auto mapIter = curParentBLP->begin(); mapIter != curParentBLP->end(); ++mapIter)
                {
                    //int radius = mapIter->first;
                    std::vector<std::vector<int>> *entity = &(mapIter->second);
                    int parentRadius = (mapIter->first);

                    int numBoundaries = (int)entity->size();//external (ancestoral boundaries)
                    for(int i=0; i<(int)numBoundaries; ++i)
                    {

                        //for(int r=-wbl; r<=wbl; ++r)
                        //for(int r=1; r<=wbl; ++r)
                        for(int r=-wbl; r<=wbl; ++r)
                        {
                            Range curRange;
                            if(r==-1)
                            {
                                //we dont need to distinguish
                                //between negative, positive and zero for external boundaries in this case,
                                //this reduces the number of regions considerably

                                r = 1;//skip other r till r=1
                                curRange.lo = (mapIter->second)[i][hopelessStart-r];
                                curRange.hi = (mapIter->second)[i][hopelessStart+r+1];
                            }
                            else
                            {
                                curRange.lo = (mapIter->second)[i][hopelessStart+r];
                                curRange.hi = (mapIter->second)[i][hopelessStart+r+1];
                            }

                            if(curRange.hi > curRange.lo)
                            {
                                //check if this is continuous to some previous
                                //range, if yes then merge 
                                bool merge= false; //needed, even if merge is switched on or off

                                //switch off merging, as it is causing problem
                                //if there is only one level
                               /* int prevRegionSize = (int)curLeaf.boundaryRange[wr][r].size();

                                for(int prevReg=0; prevReg<prevRegionSize; ++prevReg)
                                {
                                    Range prevRange = curLeaf.boundaryRange[wr][r][prevReg];
                                    if(curRange.lo == prevRange.hi)
                                    {
                                        //merge
                                        curLeaf.boundaryRange[wr][r][prevReg].hi = curRange.hi;
                                        merge = true;
                                    }
                                    else if(curRange.hi == prevRange.lo)
                                    {
                                        //merge
                                        curLeaf.boundaryRange[wr][r][prevReg].lo = curRange.lo;
                                        merge = true;
                                    }
                                }*/

                                if(!merge)
                                {
                                    int cur_wr = std::max(wr, std::abs(parentRadius)-1);
                                    curLeaf.boundaryRange[cur_wr][r].push_back(curRange);
                                }
                            }
                            else
                            {
                                printf("omiting [%d, %d]  boundary\n", curRange.lo, curRange.hi);
                            }
                        }
                    }
                }
            }
        }

        //push current boundary
        //std::vector<std::map<int,Range>> curBoundaryRange(wbl);
        for(int r=-wbl; r<=wbl; ++r)
        {
            if(r!=0)
            {
                Range curRange;
                curRange.lo = parentLeaf.lp[hopelessStart+r];
                curRange.hi = parentLeaf.lp[hopelessStart+r+1];
                //push only if it is non empty
                if(curRange.hi > curRange.lo)
                {
                    curLeaf.boundaryRange[std::abs(r)-1][r].push_back(curRange);//push direct boudary
                }
            }
        }

        printf("@@@ range = %d,%d\n", curLeaf.range.lo, curLeaf.range.hi);
        //TODO: mtxPower with boundaryRange
        mtxPower curStage(graph, highestPower, 1, cacheSize, safetyFactor, get_cache_violation_cutoff(curLeaf.stage), curLeaf.range.lo, curLeaf.range.hi, curLeaf.boundaryRange, curLeaf.nodeId, numSharedCache);
        curStage.findPartition();

        int *perm_curStage;
#ifdef RACE_PERMUTE_ON_FLY
        READ_STAGE_DATA(curLeaf, curStage, false);
#else
        READ_STAGE_DATA(curLeaf, curStage, true);
        //permute
        updatePerm(&perm, perm_curStage, totalRows, totalRows);
        delete[] perm_curStage;
#endif

       //push leaf to tree
        tree.push_back(curLeaf);
        //make this a child of parent
        int curIdx = (int)tree.size()-1;
        tree[parentIdx].children.push_back(curIdx);
        if((push_nodeId+1) < totalNodes)
        {
            tree[parentIdx].childrenNodeStart[push_nodeId+1]++;
        }
        if(!curLeaf.hid.empty())
        {
            printf("recursively calling for parent = %d\n", curIdx);
            recursivePartition(curIdx);
        }
    }

    //sum-up the buckets of childrenNodeStart
    for(int node=1; node<totalNodes; ++node)
    {
        tree[parentIdx].childrenNodeStart[node] += tree[parentIdx].childrenNodeStart[node-1];
    }
}

void mtxPowerRecursive::findPartition()
{
    MPLeaf curLeaf;

    std::vector<int> *distFromRemotePtr = getDistFromRemotePtr();

    int numChunks = (distFromRemotePtr->size() == 0) ? 1 : (distFromRemotePtr->size() - 1);
    int mainRowsBegin = (distFromRemotePtr->size() == 0) ? 0 : distFromRemotePtr->at(numChunks - 1);
    int mainRowsEnd = (distFromRemotePtr->size() == 0) ? -1 : distFromRemotePtr->at(numChunks);
    bool haveMPI = !distFromRemotePtr->empty();

    bool printCheck = true;
    if(printCheck && haveMPI){
        printf("\n-------------- distFromRemotePtr check --------------\n");
        printf("distFromRemotePtr->size() = %i\n", distFromRemotePtr->size());
        if(distFromRemotePtr->size() != 0) {
            printf("distFromRemotePtr = [");
            for(int i = 0; i < distFromRemotePtr->size(); ++i){
                if(i == (distFromRemotePtr->size() - 1)){
                    printf("%i", distFromRemotePtr->at(i));
                    break;
                }
                printf("%i, ", distFromRemotePtr->at(i));
            }
            printf("]\n");
        }
        printf("numChunks = %i\n", numChunks);
        printf("mainRowsBegin = %i\n", mainRowsBegin);
        printf("mainRowsEnd = %i\n", mainRowsEnd);
        printf("-----------------------------------------------------\n\n");
    }
    curLeaf.parent = -1;
    curLeaf.nodeId = -1;
    curLeaf.range.lo = mainRowsBegin;
    curLeaf.range.hi = (mainRowsEnd == -1) ? graph->NROW : mainRowsEnd;
    curLeaf.stage = 0;
    curLeaf.boundaryRange = {};
// TODO: include in config.h.in
// #ifdef RACE_HAVE_MPI 
    if(haveMPI){
        int wbl = workingBoundaryLength_base(highestPower);
        curLeaf.boundaryRange.resize(wbl);
        for(int r=-wbl; r<0; ++r) //only -ve since only input dependencies
        {
            Range curRange;
            curRange.lo = distFromRemotePtr->at(numChunks - 1 + r);
            curRange.hi = distFromRemotePtr->at(numChunks + r);
            if(printCheck){
                printf("----------------------------\n");
                printf("curRange.lo = %i\n", curRange.lo);
                printf("curRange.hi = %i\n", curRange.hi);
                printf("----------------------------\n");
            }
            //push only if it is non empty
            if(curRange.hi > curRange.lo)
            {
                curLeaf.boundaryRange[std::abs(r)-1][r].push_back(curRange);//push direct boudary
            }
        }
    }
// #endif
    //partition for first stage
    mtxPower curStage(graph, highestPower, numSharedCache, cacheSize, safetyFactor, get_cache_violation_cutoff(curLeaf.stage), curLeaf.range.lo, curLeaf.range.hi, curLeaf.boundaryRange, -1, -1, mtxType);
    curStage.findPartition();
    hopelessNodePtr = curStage.getHopelessNodePtr();
    int* perm_curStage;
#ifdef RACE_PERMUTE_ON_FLY
    READ_STAGE_DATA(curLeaf, curStage, false);
#else
    READ_STAGE_DATA(curLeaf, curStage, true);
    updatePerm(&perm, perm_curStage, totalRows, totalRows);
    delete[] perm_curStage;
#endif

    tree.push_back(curLeaf);

    if(mtxType == "N")
    {
        if(!curLeaf.hid.empty())
        {
            recursivePartition(0); //my curId=0
        }
    }

#ifndef RACE_PERMUTE_ON_FLY
    UPDATE_INVPERM(invPerm, perm);
#else
    int len;
    graph->getPerm(&perm, &len);
    graph->getInvPerm(&invPerm, &len);
#endif
    //totalLevel = (int)(curLeaf.lp.size()-1);
    //if empty this is final
    /*COPY_TO_PLAIN_INT_PTR(curLeaf.lp, lp);
    COPY_TO_PLAIN_INT_PTR(nodePtr_vec, nodePtr);
    COPY_TO_PLAIN_INT_PTR(curLeaf.unlockRow, unlockRow);
    COPY_TO_PLAIN_INT_PTR(curLeaf.unlockCtr, unlockCtr);
    COPY_TO_PLAIN_INT_PTR(curLeaf.dangerRow, dangerRow);
*/

    printTree();
}

void mtxPowerRecursive::printTree()
{
    for(int i=0; i<(int)tree.size(); ++i)
    {
        MPLeaf curLeaf = tree[i];
        printf("%d. Range:[%d, %d], ", i, curLeaf.range.lo, curLeaf.range.hi);
        printf("Children:[");
        for(int j=0; j<(int)curLeaf.children.size(); ++j)
        {
            printf("%d, ", curLeaf.children[j]);
        }
        printf("], ChildrenNodeStart:[");
        for(int j=0; j<(int)curLeaf.childrenNodeStart.size(); ++j)
        {
            printf("%d, ", curLeaf.childrenNodeStart[j]);
        }
        printf("], Parent: %d, node: %d, stage: %d, cache_cutoff: %d, hopelessRegions: [",
                curLeaf.parent, curLeaf.nodeId, curLeaf.stage, get_cache_violation_cutoff(curLeaf.stage));
        for(int j=0; j<(int)curLeaf.hid.size(); ++j)
        {
            printf("%d:%d, ", curLeaf.lp[curLeaf.hid[j]], curLeaf.lp[curLeaf.hid[j]+1]);
        }
        printf("], hid: [");
        for(int j=0; j<(int)curLeaf.hid.size(); ++j)
        {
            printf("%d ", curLeaf.hid[j]);
        }
        printf("], unitPtr: [");
        for(int j=0; j<(int)curLeaf.unitPtr.size(); ++j)
        {
            printf("%d ", curLeaf.unitPtr[j]);
        }
        printf("], unitNodePtr: [");
        for(int j=0; j<(int)curLeaf.unitNodePtr.size(); ++j)
        {
            printf("%d ", curLeaf.unitNodePtr[j]);
        }
        printf("], Boundaries: [");
        EXEC_BOUNDARY_STRUCTURE(curLeaf.boundaryRange,
                printf("{%d,%d} -> (%d, %d) ", _workingRadius_, _radius_, _val_.lo, _val_.hi);
                );
        printf("]\n");
    }
}

/*int mtxPowerRecursive::getTotalLevel()
{
    return totalLevel;
}
*/
/*
int mtxPowerRecursive::getTotalNodes()
{
    return numSharedCache;
}
*/
/*
int*  mtxPowerRecursive::getLevelPtrRef()
{
    return lp;
}

int*  mtxPowerRecursive::getUnlockRowRef()
{
    return unlockRow;
}


int*  mtxPowerRecursive::getDangerRowRef()
{
    return dangerRow;
}

int*  mtxPowerRecursive::getUnlockCtrRef()
{
    return unlockCtr;
}
*/
/*
std::vector<int> getNodePtr()
{
    return nodePtr;
}
*/
/*
std::vector<int> getHopelessNodePtr()
{
    return hopelessNodePtr;
}*/

void mtxPowerRecursive::getPerm(int **perm_, int *len)
{
    (*perm_) = perm;
    perm = NULL; //hand-over responsibility
    (*len) = totalRows;
}

void mtxPowerRecursive::getInvPerm(int **invPerm_, int *len)
{
    (*invPerm_) = invPerm;
    invPerm = NULL; //hand-over responsibility
    (*len) = totalRows;
}

std::vector<int>* mtxPowerRecursive::getDistFromRemotePtr()
{
    return &(graph->distFromRemotePtr);
}
