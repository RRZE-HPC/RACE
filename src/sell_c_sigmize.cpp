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

#include "sell_c_sigmize.h"

inline void sell_c_sigmize_Kernel(int start, int end, void *args)
{
    sell_c_sigmize_arg *mat = (sell_c_sigmize_arg *)(args);
    int simdWidth = mat->simdWidth;

    int C = mat->C;
    for(int row=start; row<end; ++row)
    {
        int chunk = row/C;
        int rowinchunk = row%C;
        int base_idx = mat->chunkStart[chunk] + rowinchunk;
        for(int j=mat->rowLen[row]; j<mat->chunkLenPadded[chunk]; ++j)
        {
            int idx = base_idx + j*C;
            //row%simdWidth is required; to avoid
            //false sharing within a simd lane;
            //now each 4 row gets different values
            //avoiding such false sharing
            mat->col[idx] = start + ((row-start)%simdWidth);
        }
    }
}

void sell_c_sigmize(int simdWidth, int C, int* col, int* chunkStart, int* rl, int *clp, RACE::Interface *ce)
{
    sell_c_sigmize_arg *mat = new sell_c_sigmize_arg;
    mat->simdWidth = simdWidth;
    mat->C = C;
    mat->col = col;
    mat->chunkStart = chunkStart;
    mat->rowLen = rl;
    mat->chunkLenPadded = clp;

    void* voidMat = (void*) (mat);
    int sell_c_sigma_id = ce->registerFunction(&sell_c_sigmize_Kernel, voidMat);
    ce->executeFunction(sell_c_sigma_id);
}
