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

#ifndef RACE_SELL_C_SIGMIZE
#define RACE_SELL_C_SIGMIZE

#include "interface.h"

struct sell_c_sigmize_arg
{
    int simdWidth;
    int C;
    int* col;
    int* chunkStart;
    int* rowLen;
    int* chunkLenPadded;
};


void sell_c_sigmize_Kernel(int start, int end, void *args);
void sell_c_sigmize(int simdWidth, int C, int* col, int* chunkStart, int* rl, int *clp, RACE::Interface* ce);

#endif
