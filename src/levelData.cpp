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

#include "levelData.h"

LevelData::LevelData():levelRow(NULL),levelNnz(NULL),totalLevel(0)
{
    nrow = 0;
    nnz = 0;
}

//assignment operator
LevelData& LevelData::operator=(const LevelData& other)
{
    if(this!=&other)
    {
        if(totalLevel != 0)
        {
            ERROR_PRINT("Cannot assign levelData")
        }
        totalLevel = other.totalLevel;
        levelRow = new int[totalLevel];
        levelNnz = new int[totalLevel];
        for(int i=0; i<totalLevel; ++i)
        {
            levelRow[i] = other.levelRow[i];
            levelNnz[i] = other.levelNnz[i];
        }

        nrow = other.nrow;
        nnz  = other.nnz;
    }

    return *this;
}


LevelData::~LevelData()
{
    if(levelRow) {
        delete[] levelRow;
    }

    if(levelNnz) {
        delete[] levelNnz;
    }
}
