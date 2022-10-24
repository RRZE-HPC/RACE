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

//if GAP is there use optimised version

#include "config.h"

#ifdef RACE_USE_GAP
    #ifdef RACE_USE_SOA_GRAPH
        #include "graph_SoA.cpp"
        #pragma message ("SoA graph being used")
    #else
        #include "graph_AoS.cpp"
        #pragma message ("AoS graph being used")
    #endif
#else
    #ifdef RACE_USE_SOA_GRAPH
        #pragma message ("SoA graph not supported with serial BFS. Switch on RACE_USE_GAP. Now AoS graph being used")
    #else
        #pragma message ("AoS graph being used")
    #endif
    #include "graph_AoS.cpp"
#endif

#include "graph.h"

RACE_error RACE::makeSymmetricGraph(int NROW, int NCOL, int* rowPtr, int* col, int **outRowPtr, int **outCol)
{
    //in-case of MPI need to adapt this
    if (NROW != NCOL) {
        ERROR_PRINT("Not square matrix. RACE currently considers only square matrices.");
        return RACE_ERR_MATRIX_SYMM;
    }

    int newNNZ=rowPtr[NROW];
    std::vector<std::vector<int>> newColInRow(NROW);

    volatile bool isSymmetric = true;
#pragma omp parallel for schedule(static)
    for(int r = 0; r <NROW; ++r)
    {
        for(int j = rowPtr[r]; j < rowPtr[r+1]; ++j)
        {
            int c = col[j];
            if(c != r && c >= 0 && c < NCOL)
            {
                bool hasPair = false;
                for(int k = rowPtr[c]; k < rowPtr[c+1]; ++k)
                {
                    if(col[k] == r)
                    {
                        hasPair = true;
                        break;
                    }
                }
                if(!hasPair)
                {
                    //add the missing pair
                    isSymmetric = false;
                    //do an update to the rows and col
#pragma omp critical
                    {
                        newColInRow[c].push_back(r);
                        ++newNNZ;
                    }
                }
            }
        } // loop over col in row
    } // loop over row

    if(!isSymmetric)
    {
        INFO_PRINT("Matrix does not have a symmetric pattern.");
        //make new data structure adding the missing entries
        (*outRowPtr) = new int[NROW+1];
        (*outCol) = new int[newNNZ];
        int *newRowPtr = (*outRowPtr);
        int *newCol = (*outCol);

        //update rowPtr
        newRowPtr[0] = rowPtr[0];
        //NUMA init
#pragma omp parallel for schedule(static)
        for(int r = 0; r <NROW; ++r)
        {
            newRowPtr[r+1] = 0;
        }
        for(int r = 0; r <NROW; ++r)
        {
            //new rowLen = old + extra nnz;
            int rowLen = (rowPtr[r+1]-rowPtr[r]) + newColInRow[r].size();
            newRowPtr[r+1] = newRowPtr[r]+rowLen;
        }

        if(newNNZ != newRowPtr[NROW])
        {
            ERROR_PRINT("Internal error: new NNZ count does not match last entry in new rowPtr");
            return RACE_ERR_INTERNAL;
        }
        //update col
#pragma omp parallel for schedule(static)
        for(int r = 0; r <NROW; ++r)
        {
            int j, old_j;
            for(j=newRowPtr[r], old_j=rowPtr[r]; old_j<rowPtr[r+1]; ++j, ++old_j)
            {
                newCol[j] = col[old_j];
            }
            //add extra nnzs now
            int ctr=0;
            for(; j< newRowPtr[r+1]; ++j)
            {
                newCol[j] = newColInRow[r][ctr];
                ++ctr;
            }
        }

        PERFWARNING_PRINT("However, RACE internally added %3.2fx extra non-zeros to convert the matrix to a symmetric pattern for computaiton of the permutation. There are better ways to do it, but we currently stick to this approach for simplicity.", (newRowPtr[NROW])/static_cast<double>(rowPtr[NROW]));
    }

    return RACE_SUCCESS;
}
