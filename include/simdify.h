#ifndef RACE_SIMDIFY_H
#define RACE_SIMDIFY_H

#include "print.h"
#include "interface.h"
#include "sell_c_sigmize.h"
#include <vector>
#include <algorithm>

template <typename T> bool simdifyTemplate(int simdWidth, int C, int nrows, int* col, int* chunkStart, int* rl, int *clp, T *val, RACE::Interface *ce);

//print column entry corresponding to the row in arg (for debugging purposes)
/*static void print_row(int row, int C,  int *chunkStart, int *col, int* clp)
{
    printf("col in row-%d :", row);
    int chunk = row/C;
    int rowinchunk = row%C;
    int idx = chunkStart[chunk] + rowinchunk;
    for(int j=0; j<clp[chunk]; ++j)
    {
        printf("%d ",col[idx]);
        idx+=32;
    }
    printf("\n");
}*/

/*@brief: This function re-arranges column indices within a row to enable
 * simd operations. This is applicable only for SELL-C-sigma formats.
 * After construction of the matrix call this function in order to make
 * SELL-C-sigma work for distance-2 dependent kernels.
 *
 * @param[in] simdWidth Width of simd vector
 * @param[in] C chunkLen, it has to be a multiple of simdWidth
 * @param[in] nrows number of rows in the matrix
 * @param[in/out] col column index of the matrix
 * @param[in] rl rowLen of the matrix
 * @param[in] cl  chunkLen of each chunk
 * @param[in] chunkStart padded chunkLen of each chunk
 * @param[in/out] val non zeros of the matrix
 *
*/
template <typename T> bool simdifyTemplate(int simdWidth, int C, int nrows, int* col, int* chunkStart, int* rl, int *clp, T *val, RACE::Interface *ce)
{
    if(C%simdWidth)
    {
        ERROR_PRINT("Chunk length(C) of the  matrix not a multiple of simd width, this should not have happened");
    }

    //sell-C-sigmize the kernel
    sell_c_sigmize(simdWidth, C, col, chunkStart, rl, clp, ce);

    int nChunk = static_cast<int>(nrows/(static_cast<double>(C)));
    //int remChunk = nrows%C;

    int simdInChunk = C/simdWidth;

    //int negativeNum = -1;

    int *currCol = NULL;
    bool repeat = false;
    int ctr = 0, exchangeIdx = 0;

    bool *collisionCtr = new bool[nrows];

    //std::vector<int> dummyCol;

    for(int i=0; i<nrows; ++i)
    {
        collisionCtr[i] = false;
    }

    bool wrapped = false;

    //for debugging
    //right now, it  just checks whether  every previous
    //column entries are still present; TODO
    bool checkCorrectness = false;

    //col indices that a simd register would
    //operate at a time
    int *simdCol = new int [simdWidth];

    std::vector<int> col_duplicate;


    for(int row=0; row<nrows; ++row)
    {
        int chunk = row/C;
        int rowinchunk = row%C;
        for(int j=0; j<clp[chunk]; ++j)
        {
            int col_idx = chunkStart[chunk] + rowinchunk + j*C;
            col_duplicate.push_back(col[col_idx]);
        }
    }

    for(int chunk=0; chunk<(nChunk+1); ++chunk)
    {
        for(int simdIdx=0; simdIdx<simdInChunk; ++simdIdx)
        {

#if 0
            //this technique need not be done
            //and would cause problems; since
            //we can't distinguish real and dummy
            //entries; therefore now the dummy has
            //to be distinctive in successive simdWidth
            //rows; and we directly apply the rest of the
            //method on it

            negativeNum = -1;
            dummyCol.clear();

            //assign dummy-values to distinctive values to make computations
            //easier
            row = chunk*simdInChunk*simdWidth + simdIdx*simdWidth;
            for(int i=0; ((i<simdWidth) && (row<nrows)); ++i, ++row)
            {
                //row = chunk*simdInChunk*simdWidth + simdIdx*simdWidth + i;
                for(int j=rl[row]; j<clp[chunk]; ++j)
                {
                    col_idx = chunkStart[chunk]+C*j+simdIdx*simdWidth+i;
                    dummyCol.push_back(col[col_idx]);
                    col[col_idx] = negativeNum;
                    --negativeNum;
                }
            }
#endif

            int row=0, col_idx=0;
            //re-arrange simdWidth rows here, in a fashion to enable SIMD
            for(int j=0; j<clp[chunk]; ++j)
            {
                row = chunk*simdInChunk*simdWidth + simdIdx*simdWidth;
                for(int i=0; ((i<simdWidth) && (row<nrows)); ++i, ++row)
                {
                    col_idx = chunkStart[chunk]+C*j+simdIdx*simdWidth+i;
                    row = chunk*simdInChunk*simdWidth + simdIdx*simdWidth + i;
                    currCol = &(col[col_idx]);
                    repeat = true;
                    ctr = 0;
                    exchangeIdx = j;

                    wrapped = false;
                    //switch with next entry; where no repetition is there
                    while(repeat)
                    {
                        repeat = false;
                        //check whether currCol is already there
                        for(int k=0; k<i; ++k)
                        {
                            //it's there, requires rearrangement
                            if( (*currCol) == simdCol[k] )
                            {
                               repeat = true;
                                collisionCtr[row] = true;
                                break;
                            }
                        }

                        //need to check for counterpart also
                        if(wrapped == true)
                        {
                            for(int k=0; k<simdWidth; ++k)
                            {
                                if(k!=i)
                                {
                                    //form exchange Column indices
                                    int exchangeCol = col[chunkStart[chunk]+C*exchangeIdx+simdIdx*simdWidth+k];
                                    //check if the exchange is compatible
                                    if(col[col_idx] == exchangeCol )
                                    {
                                        repeat = true;
                                        break;
                                    }
                                }
                            }
                        }

                        //printf("repeat = %d\n", repeat);
                        if(repeat)
                        {
                            if(ctr < clp[chunk])
                            {
                                if(exchangeIdx != (clp[chunk]-1))
                                {
                                    currCol += C;
                                    ++exchangeIdx;
                                    //printf("exchange idx = %d, col = %d\n", exchangeIdx, *currCol);
                                }
                                //wrap around
                                else
                                {
                                    currCol = &(col[chunkStart[chunk]+simdIdx*simdWidth+i]);
                                    exchangeIdx = 0;
                                    wrapped = true;
                                }
                                ++ctr;
                            }
                            else
                            {
                                ERROR_PRINT("Cannot simdify check whether a min. rowLen of %d is present rowLen = %d chunkLen = %d row =  %d", simdWidth, rl[row], clp[chunk], row);
                                return false;
                            }
                        }
                    }
                    simdCol[i] = (*currCol);

                    //exchange current col and nnz with exchangeIdx-th element
                    if(exchangeIdx != j)
                    {
                        int tempCol = col[chunkStart[chunk]+simdIdx*simdWidth+C*j+i];
                        col[chunkStart[chunk]+simdIdx*simdWidth+C*j+i] = col[chunkStart[chunk]+simdIdx*simdWidth+C*exchangeIdx+i];
                        col[chunkStart[chunk]+simdIdx*simdWidth+C*exchangeIdx+i] = tempCol;

                        T temp_val = val[chunkStart[chunk]+simdIdx*simdWidth+C*j+i];
                        val[chunkStart[chunk]+simdIdx*simdWidth+C*j+i] = val[chunkStart[chunk]+simdIdx*simdWidth+C*exchangeIdx+i];
                        val[chunkStart[chunk]+simdIdx*simdWidth+C*exchangeIdx+i] = temp_val;
                    }
                //    printf("%d ", col[chunkStart[chunk]+C*j+i]);
                }
               // printf("\n");
            }

#if 0
            int dummy_ctr = 0;
            //re-assign negative values
            for(int j=0; j<clp[chunk]; ++j)
            {
                row = chunk*simdInChunk*simdWidth + simdIdx*simdWidth;
                for(int i=0; ((i<simdWidth)&&(row<nrows)); ++i, ++row)
                {
                    col_idx = chunkStart[chunk]+C*j+simdIdx*simdWidth+i;
/*
                    if(row==2)
                        printf("after: %d ", col[col_idx]);
*/

                    if(col[col_idx] < 0)
                    {
                        if(val[col_idx] != 0)
                        {
                            ERROR_PRINT("ERROR in simdify");
                        }
                        //dummy_col mapped to original index
                        col[col_idx] = dummyCol[(-col[col_idx])-1];
                        dummy_ctr++;
                    }

                }
            }

            if(dummy_ctr != dummyCol.size())
            {
                ERROR_PRINT("Internal Error: mismatch in dummy columns");
            }
#endif

        }
    }


    //count collisions to provide statistics
    int totalCollision = 0;
    for(int i=0; i<nrows; ++i)
    {
        if(collisionCtr[i] == true)
            ++totalCollision;
    }

    int start_nnz = 0;
    int nnz = 0;

    if(checkCorrectness)
    {
        //check correctness
        for(int row=0; row<nrows; ++row)
        {
            int chunk = row/C;
            int rowinchunk = row%C;
            start_nnz = nnz;

            //copy
            std::vector<int> find_in_col(clp[chunk]);

            for(int j=0; j<clp[chunk]; ++j)
            {
                int col_idx = chunkStart[chunk] + rowinchunk + j*C;
                auto first = col_duplicate.begin()+start_nnz;
                auto last = col_duplicate.begin()+start_nnz+clp[chunk];
                auto it = std::find(first, last, col[col_idx]);
                if(it == last)
                {
                    ERROR_PRINT("\nNon-compatible permutations in simdify, Check row = %d could not find col = %d rowLen = %d chunkLen = %d", row, col[col_idx], rl[row], clp[chunk]);
                    printf("old col =\n");
                    for(int k=0; k<clp[chunk]; ++k)
                    {
                        printf("%d ",col_duplicate[start_nnz+k]);
                    }
                    printf("\n\nnew col =\n");
                    for(int k=0; k<clp[chunk]; ++k)
                    {
                        int idx = chunkStart[chunk] + rowinchunk + k*C;
                        printf("%d ",col[idx]);
                    }
                    printf("\n");

                    exit(0);
                }
                ++nnz;
            }
        }
    }

    INFO_PRINT("SIMDIFY: resolved %d row conflicts for simdWidth = %d", totalCollision, simdWidth);


    delete[] collisionCtr;
    delete[] simdCol;

    return true;
}


bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int*rl, int* clp, double* val);

bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int*rl, int* clp, float* val);

bool testSimdify();

#endif
