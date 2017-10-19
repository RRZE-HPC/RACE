#include "sparsemat.h"
#include "mmio.h"
#include "stdlib.h"
#include <omp.h>
#include <vector>

sparsemat::sparsemat():nrows(0), nnz(0), ce(NULL), val(NULL), rowPtr(NULL), col(NULL), nnz_symm(0), rowPtr_symm(NULL), col_symm(NULL), val_symm(NULL)
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

//necessary for GS like kernels
void sparsemat::makeDiagFirst()
{
    //check whether a new allocation is necessary
    int extra_nnz=0;
    std::vector<double>* val_with_diag = new std::vector<double>;
    std::vector<int>* col_with_diag = new std::vector<int>;

    for(int row=0; row<nrows; ++row)
    {
        bool diagHit = false;
        for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
        {
            val_with_diag->push_back(val[idx]);
            col_with_diag->push_back(col[idx]);

            if(col[idx] == row)
            {
                diagHit = true;
                break;
            }
        }
        if(!diagHit)
        {
            val_with_diag->push_back(0.0);
            col_with_diag->push_back(row);
            ++extra_nnz;
            rowPtr[row+1] = rowPtr[row+1] + extra_nnz;
        }
    }

    //allocate new matrix if necessary
    if(extra_nnz)
    {
        delete[] val;
        delete[] col;

        nnz += extra_nnz;
        val = new double[nnz];
        col = new int[nnz];

        for(int idx=0; idx<nnz; ++idx)
        {
            val[idx] = val_with_diag->at(idx);
            col[idx] = col_with_diag->at(idx);
        }
        printf("Explicit 0 in diagonal entries added\n");
    }

    delete val_with_diag;
    delete col_with_diag;

#pragma omp parallel for
    for(int row=0; row<nrows; ++row)
    {
        bool diag_hit = false;

        double* newVal = new double[rowPtr[row+1]-rowPtr[row]];
        int* newCol = new int[rowPtr[row+1]-rowPtr[row]];
        for(int idx=rowPtr[row], locIdx=0; idx<rowPtr[row+1]; ++idx, ++locIdx)
        {
            //shift all elements+1 until diag entry
            if(col[idx] == row)
            {
                newVal[0] = val[idx];
                newCol[0] = col[idx];
                diag_hit = true;
            }
            else if(!diag_hit)
            {
                newVal[locIdx+1] = val[idx];
                newCol[locIdx+1] = col[idx];
            }
            else
            {
                newVal[locIdx] = val[idx];
                newCol[locIdx] = col[idx];
            }
        }
        //assign new Val
        for(int idx = rowPtr[row], locIdx=0; idx<rowPtr[row+1]; ++idx, ++locIdx)
        {
            val[idx] = newVal[locIdx];
            col[idx] = newCol[locIdx];
        }

        delete[] newVal;
        delete[] newCol;
    }
}

//write matrix market file
bool sparsemat::writeFile(char* filename)
{
    int* row = new int[nnz];
    //create row indices
    for(int i=0; i<nrows; ++i)
    {
        for(int idx=rowPtr[i]; idx<rowPtr[i+1]; ++idx)
        {
            row[idx]=i+1;
            col[idx]+=1;
        }
    }

    mm_write_mtx_crd(filename, nrows, nrows, nnz, row, col, val, "MCRG");

    delete[] row;
}

bool sparsemat::computeSymmData()
{
    //compute only if previously not computed
    if(nnz_symm == 0)
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
    }
    return true;
}

void sparsemat::colorAndPermute(dist_t dist, int nthreads, int smt, PinMethod pinMethod)
{
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

void sparsemat::NUMA_init(bool symmPart)
{
    double* targetVal = (symmPart)?val_symm:val;
    int* targetCol = (symmPart)?col_symm:col;
    int* targetRowPtr = (symmPart)?rowPtr_symm:rowPtr;
    int targetNnz = (symmPart)?nnz_symm:nnz;

    double* NUMA_val = new double[targetNnz];
    int* NUMA_col = new int[targetNnz];
    int* NUMA_rowPtr = new int[nrows+1];

    NUMA_rowPtr[0] = targetRowPtr[0];

#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        NUMA_rowPtr[row+1] = targetRowPtr[row+1];
        for(int idx=targetRowPtr[row]; idx<targetRowPtr[row+1]; ++idx)
        {
            NUMA_val[idx] = targetVal[idx];
            NUMA_col[idx] = targetCol[idx];
        }
    }

    //now delete old arrays
    delete[] targetVal;
    delete[] targetCol;
    delete[] targetRowPtr;

    if(symmPart)
    {
        val_symm = NUMA_val;
        col_symm = NUMA_col;
        rowPtr_symm = NUMA_rowPtr;
    }
    else
    {
        val = NUMA_val;
        col = NUMA_col;
        rowPtr = NUMA_rowPtr;
    }
}
