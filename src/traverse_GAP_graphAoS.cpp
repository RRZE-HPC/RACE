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

#include "traverse.h"
#include "utility.h"
#include <set>
#include <omp.h>
#include "timing.h"

std::map<int, LevelData> RACE::Traverse::cachedData;
RACE::Traverse::Traverse(RACE::Graph *graph_, RACE::dist dist_, int rangeLo_, int rangeHi_, int parentIdx_, int numRoots_, std::vector<std::map<int, std::vector<Range>>> boundaryRange_, str::string mtxType_):graph(graph_),dist(dist_), rangeLo(rangeLo_),rangeHi(rangeHi_),parentIdx(parentIdx_), numRoots(numRoots_), graphSize(graph_->graphData.size()),distFromRoot(NULL),perm(NULL),invPerm(NULL), boundaryRange(boundaryRange_), boundary_bm(NULL), queue(graphSize), levelData(NULL), mtxType(mtxType_)
{
    if( (mtxType != "N") && ( (mtxType != "L" && mtxType != "U") ) )
    {

        ERROR_PRINT("Matrix type %s does not exist. Available options are: N, L, or U", mtxType.c_str());
        return;
    }


    colRangeLo = rangeLo;
    colRangeHi = rangeHi;

    if(rangeHi == -1)
    {
        rangeHi = graph->NROW;
    }

    totalThreads = omp_get_num_threads();

    distFromRoot = new int[graphSize];
    parent = new int[graphSize];
    for(int i=0; i<graphSize; ++i) {
        distFromRoot[i] = -1;
        parent[i] = -1;
    }

    //this is to store permutation vector for current BFS step only
    perm = new int[graph->NROW];
    invPerm = new int[graph->NROW];
    for(int i=0; i<graph->NROW; ++i) {
        perm[i] = i;
        invPerm[i] = i;
    }


    totalNodesIncBoundary=(rangeHi-rangeLo);//the main nodes
    minBoundary = rangeLo;
    maxBoundary = rangeHi;
    if(dist == RACE::POWER)
    {
        //if this is recursive part
        if((rangeLo!=0) || (rangeHi!=graph->NROW))
        {
            INIT_BOUNDARY_STRUCTURE(boundaryRange, boundaryLevelData, NULL);
            EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange, totalNodesIncBoundary += (_val_.hi - _val_.lo); minBoundary=std::min(minBoundary, _val_.lo); maxBoundary=std::max(maxBoundary, _val_.hi);)
            boundary_bm = new Bitmap((maxBoundary-minBoundary)+1);
            boundary_bm->reset();
            //printf("offset Boundary = %d\n", minBoundary);
            EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange, for(int node=_val_.lo; node<_val_.hi; ++node) { boundary_bm->set_bit(static_cast<size_t>(node-minBoundary));} );//use set_bit_atomic if parallelized
        }
    }

    //printf("Num boundary regions = %d\n", (int)boundaryRange.size());
    //for main (target) region
    levelData = new LevelData;
    UNUSED(parentIdx);
}

RACE::Traverse::~Traverse()
{
    if(distFromRoot) {
        delete[] distFromRoot;
    }

    if(parent)
    {
        delete[] parent;
    }

    if(boundary_bm)
    {
        delete boundary_bm;
    }

    if(perm) {
        delete[] perm;
    }

    if(invPerm) {
        delete[] invPerm;
    }

    if(levelData) {
        delete levelData;
    }

    if(dist == RACE::POWER)
    {
        EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryLevelData,
                if(_val_!=NULL)
                {
                    delete _val_;
                }
        );
    }
}

void RACE::Traverse::TDStep(int currLvl)
{
    int localCtr=0;
    int localColRangeLo = colRangeLo;
    int localColRangeHi = colRangeHi;
    //int curThreads = std::min(totalThreads, (int)(queue.size()/100.0)); //at least 100 work per thread
#pragma omp parallel //num_threads(curThreads)
    {
        QueueBuffer<int> lqueue(queue);
#pragma omp for reduction(+:localCtr) reduction(min:localColRangeLo) reduction(max:localColRangeHi)
        for( auto q_iter=queue.begin(); q_iter<queue.end(); q_iter++)
        {
            int u = *q_iter;
            int perm_u = u;
#ifdef RACE_PERMUTE_ON_FLY
            perm_u = graph->totalPerm[u];
#endif
            std::vector<int> *children = &(graph->graphData[perm_u].children);
            for(auto child_iter=children->begin(); child_iter!=children->end(); ++child_iter)
            {
                int child = *child_iter;
#ifdef RACE_PERMUTE_ON_FLY
                child = graph->totalInvPerm[child];
#endif

                //actually not needed in case of dist ONE, because you dont need
                //to worry about boundary. But doesn't hurt
                localColRangeLo = std::min(localColRangeLo, child);
                localColRangeHi = std::max(localColRangeHi, child);

                int curr_val = parent[child];
                if(curr_val == -1)
                {
                    if( (child>=rangeLo) && (child<rangeHi) )
                    {
                        if(compare_and_swap(parent[child], curr_val, u))
                        {
                            lqueue.push_back(child);
                            distFromRoot[child] = currLvl;
                            localCtr++;
                        }
                    }
                    else if((dist==RACE::POWER) && boundary_bm)
                    {
                        if((child >= minBoundary) && (child <= maxBoundary))
                        {
                            if( boundary_bm->get_bit(static_cast<size_t>(child-minBoundary)) )
                            {
                                if(compare_and_swap(parent[child], curr_val, u))
                                {
                                    lqueue.push_back(child);
                                    distFromRoot[child] = currLvl;
                                }
                            }
                        }
                    }
                    else if(dist==RACE::TWO)
                    {
                        //parent in range, then push to frontier
                        if( (u>=rangeLo) && (u<rangeHi) )
                        {
                            if(compare_and_swap(parent[child], curr_val, u))
                            {
                                lqueue.push_back(child);
                            }
                        }
                    }

                }
            }
        }
        lqueue.flush();
    }
    Counter::val = Counter::val+localCtr;
    colRangeLo = localColRangeLo;
    colRangeHi = localColRangeHi;
}

int Counter::val = 0;

void Counter::add()
{
    val++;
}

void Counter::reset()
{
    val = 0;
}

void RACE::Traverse::calculateDistance()
{
    if(mtxType == "N")
    {
        //traverse only if level has not been cached
        /*    if(cachedData.find(parentIdx) != cachedData.end())	
              {
              printf("Retrieving from cache\n");
              (*levelData) = cachedData[parentIdx];
              }
              else*/
        //bool marked_all = false;
        int root = rangeLo;
        colRangeLo = rangeLo;
        colRangeHi = rangeHi;

        Counter::reset();

        int currLvl = 0;
        for(int i=0; i<numRoots; ++i)
        {
            queue.push_back(root+i);
            distFromRoot[root+i] = currLvl;
            parent[root+i] = 0;
            Counter::add();
        }
        queue.slide_window();
        int prevIslandStop = rangeLo;

        while(!queue.empty())
        {
            currLvl += 1;
            TDStep(currLvl);
            queue.slide_window();
            //island detection
            {
                if( queue.empty() && (Counter::val != (rangeHi-rangeLo)) ) {
                    //now process islands
                    for(int i=prevIslandStop; i<rangeHi; ++i) {
                        if(distFromRoot[i] == -1) {
                            //Found him, mark him as distance 2 apart
                            if(dist != RACE::POWER) //do this only for coloring, else it might make levels with 0 elements for power kernel
                            {
                                currLvl += 1;
                            }
                            queue.push_back(i);
                            Counter::add();
                            queue.slide_window();
                            distFromRoot[i] = currLvl;
                            parent[i] = i;
                            colRangeLo = std::min(colRangeLo, i);
                            colRangeHi = std::max(colRangeHi, i);
                            prevIslandStop=i;
                            break;
                        }
                    }
                    if(queue.empty())
                    {
                        printf("counter_val = %d, range_val = %d\n", Counter::val, (rangeHi-rangeLo));
                        ERROR_PRINT("Some nodes went missing");
                        exit(-1);
                    }
                }
            }
        }

        levelData->totalLevel = currLvl;
    }
    else
    {
        if(mtxType == "L")
        {
            //forward
            for(int i=0; i<graph->NROW; ++i)
            {
                distFromRoot[i] = i;
            }
        }
        else if(mtxType == "U")
        {
            //backward
            for(int i=0; i<graph->NROW; ++i)
            {
                distFromRoot[i] = (graph->NROW-1)-i;
            }
        }
        levelData->totalLevel = graph->NROW;
    }


    printf("Total Level = %d\n",levelData->totalLevel);
    createLevelData();
    printf("created Level Data\n");
    permuteGraph();
    printf("permuted graph\n");

}

RACE_error RACE::Traverse::findLevelData(int lower_nrows, int upper_nrows, int totalLevel, LevelData* curLevelData)
{
    curLevelData->totalLevel = totalLevel;
    int* levelRow_ = new int[totalLevel];
    int* levelNnz_ = new int[totalLevel];
    for(int i=0; i<totalLevel; ++i) {
        levelRow_[i] = 0;
        levelNnz_[i] = 0;
    }

    for(int i=lower_nrows; i<upper_nrows; ++i)
    {
        int curr_dist = distFromRoot[i];
        if(curr_dist == -1)
        {
            ERROR_PRINT("There are orphan nodes; I thought this wouldn't happen");
            return RACE_ERR_GRAPH_TRAVERSAL;

        }

        levelRow_[curr_dist]+=1;
#ifdef RACE_PERMUTE_ON_FLY
        levelNnz_[curr_dist] += graph->graphData[graph->totalPerm[i]].children.size();
#else
        levelNnz_[curr_dist] += graph->graphData[i].children.size();
#endif
        //levelNnz_[curr_dist] += graph->graphData[i].upperNnz;

    }

    curLevelData->levelRow = levelRow_;
    curLevelData->levelNnz = levelNnz_;

    curLevelData->nrow=0;
    curLevelData->nnz=0;

    for(int i=0; i<totalLevel; ++i)
    {
        curLevelData->nrow += levelRow_[i];
        curLevelData->nnz += levelNnz_[i];
    }
    return RACE_SUCCESS;
}

RACE_error RACE::Traverse::createLevelData()
{
    RACE_error err_flag = RACE_SUCCESS;
    bool untouchedBoundaryNodes = false;
    //handle bondary levels in case its power calculation
    if(dist == RACE::POWER)
    {
        //set distFromRoot at active boundaries to totalLevel+1, so if it is not
        //touched from main region, it becomes the last level
        EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange,
                for(int node=_val_.lo; node<_val_.hi; ++node)
                {
                    if(distFromRoot[node] == -1)
                    {
                        distFromRoot[node] = levelData->totalLevel;
                        untouchedBoundaryNodes = true;
                    }
                }
            );

        //increase total levels if there are untouched boundary nodes
        if(untouchedBoundaryNodes)
        {
            levelData->totalLevel += 1;
        }
        EXEC_BOUNDARY_STRUCTURE(boundaryLevelData,
                _val_ = new LevelData;
                err_flag = RACE_SUCCESS;
                Range curRange = boundaryRange[_workingRadius_][_radius_][_region_];
                err_flag = findLevelData(curRange.lo, curRange.hi, levelData->totalLevel, _val_);
                if(err_flag != RACE_SUCCESS)
                {
                    ERROR_PRINT("Something went wrong in levelData calculation for boundaries");
                    return err_flag;
                }
                boundaryLevelData[_workingRadius_][_radius_][_region_] = _val_;
            );
    }

    //levelData for main body (region)
    err_flag = findLevelData(rangeLo, rangeHi, levelData->totalLevel, levelData);
    if(err_flag != RACE_SUCCESS)
    {
        ERROR_PRINT("Something went wrong in levelData calculation");
        return err_flag;
    }


    //cache this data for later use
    //Don't cache this won't happen in 
    //the current strategy
    //cachedData[parentIdx] = (*levelData);
    return err_flag;
}

void RACE::Traverse::permuteGraph()
{
    int numRegions = 1;
    std::vector<int> regionRange;
    regionRange.push_back(rangeLo);
    regionRange.push_back(rangeHi);

#ifndef RACE_PERMUTE_ON_FLY
    RACE::Graph permutedGraph(*(graph));
#endif

    if(dist==RACE::POWER)
    {
        EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange, regionRange.push_back(_val_.lo); regionRange.push_back(_val_.hi); numRegions++;);
    }

    for(int regionIdx=0; regionIdx<numRegions; ++regionIdx)
    {
        int targetRangeLo = regionRange[2*regionIdx];
        int targetRangeHi = regionRange[2*regionIdx+1];

        //create permutation vector First
        sortPerm(distFromRoot, perm, targetRangeLo, targetRangeHi);
       //create invPerm
#pragma omp parallel for schedule(static)
        for(int i=targetRangeLo; i<targetRangeHi; ++i) {
            invPerm[perm[i]] = i;
        }

#ifndef RACE_PERMUTE_ON_FLY
        //Permute rows
#pragma omp parallel for schedule(static)
        for(int i=targetRangeLo; i<targetRangeHi; ++i)
        {
            permutedGraph.graphData[i].children = graph->graphData[perm[i]].children;
        }
#endif

    }
#ifdef RACE_PERMUTE_ON_FLY
        updatePerm(&graph->totalPerm, perm, graph->NROW, graph->NROW);

#pragma omp parallel for schedule(static)
        for(int i=0; i<graph->NROW; ++i) {
            graph->totalInvPerm[graph->totalPerm[i]] = i;
        }

#endif

    if(colRangeHi != rangeHi)
    {
        colRangeHi += 1;//because we need to go 1+
    }
    if((colRangeLo > rangeLo) || (colRangeHi < rangeHi))
    {
        ERROR_PRINT("Error col range less than row range");
    }
    printf("Range = [%d,%d], colRange = [%d, %d]\n", rangeLo, rangeHi, colRangeLo, colRangeHi);

#ifndef RACE_PERMUTE_ON_FLY
    //Permute columns
    //TODO only neighbors
    //for(int i=0/*rangeLo*/; i<graph->NROW/*rangeHi*/; ++i)
#pragma omp parallel for schedule(static)
    for(int i=colRangeLo; i<colRangeHi; ++i)
    {
        std::vector<int> *children = &(permutedGraph.graphData[i].children);
        for(int j=0; j<children->size(); ++j)
        {
            int child = children->at(j);
#if 1
            bool inNodesIncBoundaries = false;
            if( (child>=rangeLo) && (child<rangeHi) ) {
                inNodesIncBoundaries = true;
            }
            if((RACE::POWER) && ((!inNodesIncBoundaries) && boundary_bm))
            {

                if((child >= minBoundary) && (child <= maxBoundary))
                {
                    if(boundary_bm->get_bit(static_cast<size_t>(child-minBoundary)))
                    {
                        inNodesIncBoundaries = true;
                    }
                }
                /*
                EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange,
                        if((child>=_val_.lo) && (child<_val_.hi))
                        {
                            inNodesIncBoundaries = true;
                            _numBoundaries_ = -1; //break from regions
                            _wbl_ = -1; //break from working radius
                            _mapIter_ = boundaryRange[_workingRadius_].end();//break from radius
                            --_mapIter_;
                            break;
                        }
                    );
                    */
            }
#endif

            if(inNodesIncBoundaries)
            {
                children->at(j) = invPerm[child];
            }
        }
    }

    permutedGraph.swap(*(graph));
#endif
}

//Getter functions
void RACE::Traverse::getPerm(int **perm_, int *len)
{
#ifndef RACE_PERMUTE_ON_FLY
    (*perm_) = perm;
    perm = NULL;
#else
    WARNING_PRINT("Use getPerm from graph directly when RACE_PERMUTE_ON_FLY is active");
    (*perm_) = graph->totalPerm;
#endif
    (*len) = graph->NROW;
}

void RACE::Traverse::getInvPerm(int **invPerm_, int *len)
{
#ifndef RACE_PERMUTE_ON_FLY
    (*invPerm_) = invPerm;
    invPerm = NULL;
#else
    WARNING_PRINT("Use getPerm from graph directly when RACE_PERMUTE_ON_FLY is active");
    (*invPerm_) = graph->totalInvPerm;
#endif
    (*len) = graph->NROW;
}

LevelData* RACE::Traverse::getLevelData()
{
    LevelData* levelData_ = levelData;
    levelData = NULL;
    return levelData_;
}

std::vector<std::map<int, std::vector<LevelData*>>> RACE::Traverse::getBoundaryLevelData()
{
    std::vector<std::map<int, std::vector<LevelData*>>> retBoundaryLevelData = boundaryLevelData;
    EXEC_BOUNDARY_STRUCTURE(boundaryLevelData,
            _val_ = NULL;
            boundaryLevelData[_workingRadius_][_radius_][_region_] = _val_;
            );
    return retBoundaryLevelData;
}
