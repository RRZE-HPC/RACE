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

#ifndef RACE_MTX_POWER_H
#define RACE_MTX_POWER_H

#include "print.h"
#include "error.h"
#include <vector>
#include "graph.h"
#include "traverse.h"
#include "levelData.h"
#include "string.h"

//this is per stage
class mtxPower
{
    RACE::Graph* graph;
    std::vector<int> nodePtr;
    std::vector<int> unlockRow;//for p2p sync
    std::vector<int> dangerRow;// for p2p sync
    std::vector<int> unlockCtr;
    int* cacheLevelGroup;
    std::vector<int> hopelessRegions; //stores level index at hopeless boundary
    std::vector<int> hopelessNodePtr;
    std::vector<std::vector<int>> hopelessRegionPositiveBoundary; //stores rows in new boundary
    std::vector<std::vector<int>> hopelessRegionNegativeBoundary; //stores rows in new boundary

    //Details of main (target) region
    int startRow, endRow;
    LevelData* levelData;
    std::vector<int> levelPtr;

    //Details of boundary of parent
    std::vector<std::map<int, std::vector<Range>>> boundaryRange;
    std::vector<std::map<int, std::vector<LevelData*>>> boundaryLevelData;
    std::vector<std::map<int, std::vector<std::vector<int>>>> boundaryLevelPtr;
    std::vector<std::map<int, std::vector<std::vector<int>>>> boundaryUnlockRow;
    std::vector<std::map<int, std::vector<std::vector<int>>>> boundaryDangerRow;
    std::vector<std::map<int, std::vector<std::vector<int>>>> boundaryUnlockCtr;

    RACE::Traverse* traverser;
    int totalLevel;
    int highestPower;
    int numSharedCache;
    double cacheSize;
    double safetyFactor;
    int cache_violation_cutoff;
    int nodeId;
    int numRootNodes;

    void identifyHopelessRegions(std::vector<int> cacheViolatedLevel);
    void getHopelessStartEnd(int count, int *start, int *end);
    void getHopelessStartEnd(int count, int *start, int *end, std::vector<int> _hopelessRegions_);
    double getElemUpperLimit(int level);
    int workingBoundaryLength();
    std::vector<int> findLevelPtr(int startNode, LevelData* curLevelData);

    public:
    //mtxType can be:
    //N: Normal
    //L: Lower traiangle
    //U: Upper triangle
    mtxPower(RACE::Graph* graph_, int highestPower_, int numSharedCache, double cacheSize_, double safetyFactor_, int cache_violation_cutoff_, int startRow_, int endRow_, std::vector<std::map<int, std::vector<Range>>> boundaryRange={}, int nodeId_=-1, int numRootNodes_=-1, std::string mtxType="N");
    ~mtxPower();
    void findPartition();
    void splitSharedCacheDomain();
    std::vector<int> findMacroLevelPtr(int* zones);
    void consolidatePartition();
    void getStatNUMA();
    void findUnlockCtr();
    void createLevelPtr();
    double getBytePerNNZ();
    void powerRun(int power, int *rowPtr, int *col, double *A, double *x);
    void getPerm(int **perm, int *len);
    void getInvPerm(int **invPerm, int *len);
    LevelData* getLevelDataRef();
    int getTotalLevel();
    int getTotalNodes();
    std::vector<int> getLevelPtr();
    std::vector<std::map<int, std::vector<std::vector<int>>>> getBoundaryLevelPtr();
    std::vector<std::map<int, std::vector<std::vector<int>>>> getBoundaryUnlockRow();
    std::vector<std::map<int, std::vector<std::vector<int>>>> getBoundaryDangerRow();
    std::vector<int> getNodePtr();
    std::vector<int> getUnlockRow();
    std::vector<int> getDangerRow();
    std::vector<int> getUnlockCtr();
    std::vector<int> getHopelessRegions();
    std::vector<int> getHopelessNodePtr();
    std::vector<std::vector<int>> getHopelessNegativeBoundaries();
    std::vector<std::vector<int>> getHopelessPositiveBoundaries();
};


#endif
