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

//file no more required added to interface
#include "simdify.h"

//generate required templates


bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int*rl, int* clp, double* val)
{
    //return testSimdify();
    return simdifyTemplate<double> (simdWidth, C, nrows, col, chunkStart, rl, clp, val);
}

bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int*rl, int* clp, float* val)
{
    return simdifyTemplate<float> (simdWidth, C, nrows, col, chunkStart, rl, clp, val);
}

bool testSimdify()
{
    int C =  4, simdWidth = 4;
    int nrows = 4;

    int *col = new int[16];
    double *val = new double[16];


    for(int i=0; i<nrows; ++i)
    {
        for(int j=0; j<nrows; ++j)
        {
            col[i*nrows+j] = i;
            val[i*nrows+j] = 1;
        }
    }

    int *chunkStart = new int[2];

    chunkStart[0] = 0;
    chunkStart[1] = 17;

    int *rl = new int[nrows];

    for(int i=0; i<nrows; ++i)
    {
        rl[i] = 4;
    }

    int *clp = new int[1];
    clp[0] = 4;

    return simdifyTemplate<double>(simdWidth, C, nrows, col, chunkStart, rl, clp, val);
}
