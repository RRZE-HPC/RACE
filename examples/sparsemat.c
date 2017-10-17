#include "sparsemat.h"
#include "mmio.h"
#include "stdlib.h"
#include <omp.h>

sparsemat::sparsemat():nrows(0), nnz(0), ce(NULL), val(NULL), rowPtr(NULL), col(NULL), val_symm(NULL), rowPtr_symm(NULL), col_symm(NULL)
{
}

sparsemat::~sparsemat()
{
    if(val)
        delete[] val;

    if(rowPtr)
        delete[] rowPtr;

    if(col)
        delete[] col;

    if(val_symm)
        delete[] val_symm;

    if(rowPtr_symm)
        delete[] rowPtr_symm;

    if(col_symm)
        delete[] col_symm;

    if(ce)
        delete ce;

}

bool sparsemat::readFile(char* filename)
{
    int *row, ncols;

    if(mm_read_unsymmetric_sparse(filename, &nrows, &ncols, &nnz, &val, &row, &col) < 0)
    {
        printf("Error in file reading\n");
        return false;
    }

    if(nrows != ncols)
    {
        printf("Currently only Symmetric matrices are supported\n");
        return false;
    }

    rowPtr = new int[nrows+1];

    int *nnzPerRow = new int[nrows];
    for(int i=0; i<nrows; ++i)
    {
        nnzPerRow[i] = 0;
    }

    //count nnz per row
    for(int i=0; i<nnz; ++i)
    {
        ++nnzPerRow[row[i]];
    }

    rowPtr[0] = 0;
    for(int i=0; i<nrows; ++i)
    {
        rowPtr[i+1] = rowPtr[i]+nnzPerRow[i];
    }

    if(rowPtr[nrows] != nnz)
    {
        printf("Error in reading matrix\n");
        return false;
    }

    return true;
}

bool sparsemat::computeSymmData()
{
    /* Here we compute symmetric data of matrix
     * which is used if necessary; upper symmetric
     * portion is stored*/

    nnz_symm = 0;
    //count non-zeros in upper-symm
    for(int i=0; i<nrows; ++i) {
        for(int j=rowPtr[i]; j<rowPtr[i+1]; ++j) {
            if(col[j]>=i) {
                ++nnz_symm;
            }
        }
    }

    rowPtr_symm = new int[nrows];
    col_symm = new int[nnz_symm];
    val_symm = new double[nnz_symm];

    rowPtr_symm[0] = 0;

    int ctr=0;
    for(int i=0; i<nrows; ++i) {
        for(int j=rowPtr[i]; j<rowPtr[i+1]; ++j) {
            if(col[j]>=i) {
                val_symm[ctr] = val[j];
                col_symm[ctr] = col[j];
                ++ctr;
            }
        }
        rowPtr_symm[i+1] = ctr;
    }

    return true;
}

void sparsemat::colorAndPermute(dist_t dist, int nthreads, int smt, PinMethod pinMethod)
{
    printf("nthreads = %d\n", nthreads);
    ce = new RACEInterface(nrows, nthreads, dist, rowPtr, col, smt, pinMethod, NULL, NULL);
    ce->RACEColor();

    int *perm, *invPerm, permLen;

    ce->getPerm(&perm, &permLen);
    ce->getInvPerm(&invPerm, &permLen);

    //now permute the matrix according to the permutation vector
    permute(perm, invPerm);

    free(perm);
    free(invPerm);
}

//symmetrically permute
void sparsemat::permute(int *perm, int*  invPerm)
{
    double* newVal = new double[nnz];
    int* newRowPtr = new int[nrows+1];
    int* newCol = new int[nnz];

    newRowPtr[0] = 0;

    int permIdx=0;
    for(int row=0; row<nrows; ++row)
    {
        //row permutation
        int permRow = perm[row];
        for(int idx=rowPtr[permRow]; idx<rowPtr[permRow+1]; ++idx)
        {
            //permute column-wise also
            newVal[permIdx] = val[idx];
            newCol[permIdx] = invPerm[col[idx]];

            //TODO: if needed sort here
            ++permIdx;
        }
        newRowPtr[row+1] = permIdx;
    }

    //free old permutations
    delete[] val;
    delete[] rowPtr;
    delete[] col;

    val = newVal;
    rowPtr = newRowPtr;
    col = newCol;
}
