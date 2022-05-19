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

#ifndef RACE_TYPE_H
#define RACE_TYPE_H

#include <vector>

namespace RACE {

    enum dist{
        ONE=1,
        TWO=2,
        POWER=3
    };

    enum d2Method{
        TWO_BLOCK,
        THREE_BLOCK /*TODO: To be done*/
    };

    enum LBTarget{
        ROW,
        NNZ
    };

    enum PinMethod{
        FILL,
        SCATTER
    };

    /*Load balncing mode, based on efficiency or just minimise effRow*/
    enum LBMode{
        MIN,
        EFFICIENCY
    };

}

struct Range
{
    int lo;
    int hi;

    Range();
};

inline Range::Range():lo(-1),hi(-1)
{
}

struct ZoneLeaf;

typedef std::vector<ZoneLeaf> tree_t;


#endif
