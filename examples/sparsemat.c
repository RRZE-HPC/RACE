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

    MM_typecode matcode;

    FILE *f;

    if ((f = fopen(filename, "r")) == NULL)
        return -1;


    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", filename);
        return -1;
    }

    fclose(f);

    bool compatible_flag = (mm_is_sparse(matcode) && mm_is_real(matcode)) && (mm_is_symmetric(matcode) || mm_is_general(matcode));

    bool symm_flag = mm_is_symmetric(matcode);

    if(!compatible_flag || symm_flag)
    {
        printf("The matrix market file provided is not supported\n");
        exit(0);
    }

    int ncols;
    int *row;
    int *col_unsorted;
    double *val_unsorted;

    if(mm_read_unsymmetric_sparse(filename, &nrows, &ncols, &nnz, &val_unsorted, &row, &col_unsorted) < 0)
    {
        printf("Error in file reading\n");
        return false;
    }
    if(nrows != ncols)
    {
        printf("Currently only Symmetric matrices are supported\n");
        return false;
    }

    //permute the col and val according to row
    int* perm = new int[nnz];
    for(int idx=0; idx<nnz; ++idx)
    {
        perm[idx] = idx;
    }

    sort_perm(row, perm, nnz);

    col = new int[nnz];
    val = new double[nnz];

    for(int idx=0; idx<nnz; ++idx)
    {
        col[idx] = col_unsorted[perm[idx]];
        val[idx] = val_unsorted[perm[idx]];
    }

    delete[] col_unsorted;
    delete[] val_unsorted;

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

    free(row);

    //writeFile("beforePerm.mtx");
    return true;
}

//necessary for GS like kernels
void sparsemat::makeDiagFirst()
{
    //check whether a new allocation is necessary
    int extra_nnz=0;
    std::vector<double>* val_with_diag = new std::vector<double>();
    std::vector<int>* col_with_diag = new std::vector<int>();
    std::vector<int>* rowPtr_with_diag = new std::vector<int>(rowPtr, rowPtr+nrows+1);

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
            }
        }
        if(!diagHit)
        {
            val_with_diag->push_back(0.0);
            col_with_diag->push_back(row);
            ++extra_nnz;
        }
        rowPtr_with_diag->at(row+1) = rowPtr_with_diag->at(row+1) + extra_nnz;
    }

    //allocate new matrix if necessary
    if(extra_nnz)
    {
        delete[] val;
        delete[] col;
        delete[] rowPtr;

        nnz += extra_nnz;
        val = new double[nnz];
        col = new int[nnz];
        rowPtr = new int[nrows+1];

        rowPtr[0] = rowPtr_with_diag->at(0);
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            rowPtr[row+1] = rowPtr_with_diag->at(row+1);
            for(int idx=rowPtr_with_diag->at(row); idx<rowPtr_with_diag->at(row+1); ++idx)
            {
                val[idx] = val_with_diag->at(idx);
                col[idx] = col_with_diag->at(idx);
            }
        }
        printf("Explicit 0 in diagonal entries added\n");
    }

    delete val_with_diag;
    delete col_with_diag;
    delete rowPtr_with_diag;

#pragma omp parallel for schedule(static)
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
    int* row_1_based = new int[nnz];
    int* col_1_based = new int[nnz];

    //create row indices
    for(int i=0; i<nrows; ++i)
    {
        for(int idx=rowPtr[i]; idx<rowPtr[i+1]; ++idx)
        {
            row_1_based[idx]=i+1;
            col_1_based[idx]=col[idx]+1;
        }
    }

    mm_write_mtx_crd(filename, nrows, nrows, nnz, row_1_based, col_1_based, val, "MCRG");

    delete[] row_1_based;
    delete[] col_1_based;
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
        rowPtr_symm = new int[nrows+1];
        rowPtr_symm[0] = 0;

        //NUMA init
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            rowPtr_symm[row+1] = 0;
        }

        //count non-zeros in upper-symm
        for(int row=0; row<nrows; ++row) {
            for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx) {
                if(col[idx]>=row) {
                    ++nnz_symm;
                }
                rowPtr_symm[row+1] = nnz_symm;
            }
        }

        col_symm = new int[nnz_symm];
        val_symm = new double[nnz_symm];

        //With NUMA init
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row) {
            int idx_symm = rowPtr_symm[row];
            for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx) {
                if(col[idx]>=row) {
                    val_symm[idx_symm] = val[idx];
                    col_symm[idx_symm] = col[idx];
                    ++idx_symm;
                }
            }
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

    //pin omp threads as in RACE for proper NUMA
    pinOMP(nthreads);

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

    //NUMA init
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        newRowPtr[row+1] = 0;
    }

    //first find newRowPtr; therefore wwe can do proper NUMA init
    int permIdx=0;
    for(int row=0; row<nrows; ++row)
    {
        //row permutation
        int permRow = perm[row];
        for(int idx=rowPtr[permRow]; idx<rowPtr[permRow+1]; ++idx)
        {
             ++permIdx;
        }
        newRowPtr[row+1] = permIdx;
    }


    //with NUMA init
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        //row permutation
        int permRow = perm[row];
        for(int permIdx=newRowPtr[row],idx=rowPtr[permRow]; permIdx<newRowPtr[row+1]; ++idx,++permIdx)
        {
            //permute column-wise also
            newVal[permIdx] = val[idx];
            newCol[permIdx] = invPerm[col[idx]];
        }
    }

    //free old permutations
    delete[] val;
    delete[] rowPtr;
    delete[] col;

    val = newVal;
    rowPtr = newRowPtr;
    col = newCol;
}

//here openMP threads are pinned according to
//RACE pinning, which is necessary for NUMA init
void sparsemat::pinOMP(int nthreads)
{
    omp_set_dynamic(0);    //  Explicitly disable dynamic teams
    omp_set_num_threads(nthreads);

#pragma omp parallel
    {
        int pinOrder = omp_get_thread_num();
        ce->pinThread(pinOrder);
    }
}


