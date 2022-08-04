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

#ifndef RACE_TRAVERSE_GAP_H
#define RACE_TRAVERSE_GAP_H
#include <vector>
#include <map>
#include "type.h"
#include "error.h"
#include "graph.h"
#include "levelData.h"
#include "macros.h"
#include "bitmap.h"
#include "sliding_queue.h"
#include <string>

//class RACE::Graph;


/**
 * @brief RACE namespace.
 */
namespace RACE
{
    class Traverse;
}

class RACE::Traverse{
    private:
        RACE::Graph *graph;
        static std::map<int, LevelData> cachedData;
        RACE::dist dist;
        int rangeLo;
        int rangeHi;
        int colRangeLo;//the min row where col in rangeLo:rangeHi are present
        int colRangeHi;//the max row where col in rangeLo:rangeHi are present
        int parentIdx;
        int numRoots;
        //Size without pure diagonal elements
        int graphSize;

        int* distFromRoot;
        int* parent;

        int* perm;
        int* invPerm;

        std::vector<std::map<int, std::vector<Range>>> boundaryRange;
        int totalNodesIncBoundary;
        int minBoundary; //min node of all boundaries
        int maxBoundary; //max node of all boundaries
        Bitmap* boundary_bm;

        SlidingQueue<int> queue;
        //Level Details
        LevelData* levelData;
        std::vector<std::map<int, std::vector<LevelData*>>> boundaryLevelData;
        int totalThreads;

        std::string mtxType;

        //std::vector<int> markChildren(int currChild, int currLvl);
        void TDStep(int curLvl);
        RACE_error findLevelData(int lower_nrows, int upper_nrows, int totalLevel, LevelData* curLevelData);
        RACE_error createLevelData();
        void permuteGraph();
    public:
        //constructor
        Traverse(RACE::Graph *graph_, RACE::dist dist, int rangeLo_=0, int rangeHi_=-1, int parentIdx=0, int numRoots=1, std::vector<std::map<int, std::vector<Range>>> boundaryRange_={}, std::string mtxType_="N");
        ~Traverse();
        void calculateDistance();

        //deletion of array's after get calling get fns is user's responsibility
        void getPerm(int  **perm_, int *len);
        void getInvPerm(int **invPerm_, int *len);
        LevelData* getLevelData();
        std::vector<std::map<int, std::vector<LevelData*>>> getBoundaryLevelData();
};

struct Counter{
    static int val;
    static void add();
    static void reset();
};

#endif
