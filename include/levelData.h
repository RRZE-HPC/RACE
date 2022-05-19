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

#ifndef RACE_LEVEL_DATA_H
#define RACE_LEVEL_DATA_H

#include "error.h"

struct LevelData{
    int *levelRow;
    int *levelNnz;
    int totalLevel;
    int nrow;
    int nnz;

    LevelData();
    LevelData& operator=(const LevelData& other);
    ~LevelData();
};

/*LevelData::LevelData():levelRow(NULL),levelNnz(NULL),totalLevel(0)
{

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
*/
#endif
