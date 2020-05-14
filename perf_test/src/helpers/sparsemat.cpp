#include "sparsemat.h"
#include "mmio.h"
#include "stdlib.h"
#include <omp.h>
#include <vector>
#include <sys/mman.h>
#include "config_eg.h"

#ifdef RACE_USE_SPMP
    #include "SpMP/CSR.hpp"
    #include "SpMP/reordering/BFSBipartite.hpp"
#endif


sparsemat::sparsemat():nrows(0), nnz(0), ce(NULL), val(NULL), rowPtr(NULL), col(NULL), nnz_symm(0), rowPtr_symm(NULL), col_symm(NULL), val_symm(NULL), rcmPerm(NULL), rcmInvPerm(NULL)
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

    if(rcmPerm)
        delete[] rcmPerm;

    if(rcmInvPerm)
        delete[] rcmInvPerm;

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

    bool compatible_flag = (mm_is_sparse(matcode) && (mm_is_real(matcode)||mm_is_pattern(matcode))) && (mm_is_symmetric(matcode) || mm_is_general(matcode));
    bool symm_flag = mm_is_symmetric(matcode);

    if(!compatible_flag)
    {
        printf("The matrix market file provided is not supported.\n Reason :\n");
        if(!mm_is_sparse(matcode))
        {
            printf(" * matrix has to be sparse\n");
        }

        if(!mm_is_real(matcode) && !(mm_is_pattern(matcode)))
        {
            printf(" * matrix has to be real or pattern\n");
        }

        if(!mm_is_symmetric(matcode) && !mm_is_general(matcode))
        {
            printf(" * matrix has to be either general or symmetric\n");
        }

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

    //If matrix market file is symmetric; create a general one out of it
    if(symm_flag)
    {
        printf("Creating a general matrix out of a symmetric one\n");

        int ctr = 0;

        //this is needed since diagonals might be missing in some cases
        for(int idx=0; idx<nnz; ++idx)
        {
            ++ctr;
            if(row[idx]!=col_unsorted[idx])
            {
                ++ctr;
            }
        }

        int new_nnz = ctr;

        int *row_general = new int[new_nnz];
        int *col_general = new int[new_nnz];
        double *val_general = new double[new_nnz];

        int idx_gen=0;

        for(int idx=0; idx<nnz; ++idx)
        {
            row_general[idx_gen] = row[idx];
            col_general[idx_gen] = col_unsorted[idx];
            val_general[idx_gen] = val_unsorted[idx];
            ++idx_gen;

            if(row[idx] != col_unsorted[idx])
            {
                row_general[idx_gen] = col_unsorted[idx];
                col_general[idx_gen] = row[idx];
                val_general[idx_gen] = val_unsorted[idx];
                ++idx_gen;
            }
        }

        free(row);
        free(col_unsorted);
        free(val_unsorted);

        nnz = new_nnz;

        //assign right pointers for further proccesing
        row = row_general;
        col_unsorted = col_general;
        val_unsorted = val_general;
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

    delete[] row;

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

    return true;
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

void sparsemat::doRCM()
{
#ifdef RACE_USE_SPMP
    int orig_threads = 1;
    printf("Doing RCM permutation\n");
#pragma omp parallel
    {
        orig_threads = omp_get_num_threads();
    }
    omp_set_num_threads(1);

    SpMP::CSR *csr = NULL;
    csr = new SpMP::CSR(nrows, nrows, rowPtr, col, val);
    rcmPerm = new int[nrows];
    rcmInvPerm = new int[nrows];
    if(csr->isSymmetric(false,false))
    {
        csr->getRCMPermutation(rcmInvPerm, rcmPerm);
    }
    else
    {
        printf("Matrix not symmetric RCM cannot be done\n");
    }
    omp_set_num_threads(orig_threads);
    delete csr;
#endif
}

int sparsemat::prepareForPower(int highestPower, int numSharedCache, double cacheSize, int nthreads, int smt, PinMethod pinMethod)
{
    ce = new Interface(nrows, nthreads, RACE::POWER, rowPtr, col, smt, pinMethod, rcmPerm, rcmInvPerm);
    ce->RACEColor(highestPower, numSharedCache, cacheSize);

    int *perm, *invPerm, permLen;
    ce->getPerm(&perm, &permLen);
    ce->getInvPerm(&invPerm, &permLen);
    permute(perm, invPerm , true);

    return 1;
}

int sparsemat::colorAndPermute(dist distance, int nthreads, int smt, PinMethod pinMethod)
{
    ce = new Interface(nrows, nthreads, distance, rowPtr, col, smt, pinMethod, rcmPerm, rcmInvPerm);
    RACE_error ret = ce->RACEColor();

    if(ret != RACE_SUCCESS)
    {
        printf("Pinning failure\n");
        return 0;
    }

    int *perm, *invPerm, permLen;

    ce->getPerm(&perm, &permLen);
    ce->getInvPerm(&invPerm, &permLen);

    //now permute the matrix according to the permutation vector
    permute(perm, invPerm);
    delete[] perm;
    delete[] invPerm;

    //pin omp threads as in RACE for proper NUMA
    pinOMP(nthreads);

    return 1;
}

double sparsemat::colorEff()
{
    return ce->getEfficiency();
}

int sparsemat::maxStageDepth()
{
    return ce->getMaxStageDepth();
}

void sparsemat::numaInit(bool RACEalloc)
{
    permute(NULL, NULL, RACEalloc);
}

//symmetrically permute
void sparsemat::permute(int *perm, int*  invPerm, bool RACEalloc)
{
    printf("NUMA allocing\n");

    double* newVal = new double[nnz];
    int* newCol = new int[nnz];
    int* newRowPtr = new int[nrows+1];

/*
    double *newVal = (double*) malloc(sizeof(double)*nnz);
    int *newCol = (int*) malloc(sizeof(int)*nnz);
    int *newRowPtr = (int*) malloc(sizeof(int)*(nrows+1));
*/

    newRowPtr[0] = 0;

    if(!RACEalloc)
    {
        //NUMA init
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            newRowPtr[row+1] = 0;
        }
    }
    else
    {
        ce->numaInitRowPtr(newRowPtr);
    }

    if(perm != NULL)
    {
        //first find newRowPtr; therefore we can do proper NUMA init
        int permIdx=0;
        printf("nrows = %d\n", nrows);
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
    }
    else
    {
        for(int row=0; row<nrows+1; ++row)
        {
            newRowPtr[row] = rowPtr[row];
        }
    }


    if(RACEalloc)
    {
        ce->numaInitMtxVec(newRowPtr, newCol, newVal, NULL);
    }

    if(perm != NULL)
    {
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
    }
    else
    {
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            for(int idx=newRowPtr[row]; idx<newRowPtr[row+1]; ++idx)
            {
                newVal[idx] = val[idx];
                newCol[idx] = col[idx];
            }
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
    int availableThreads = ce->getNumThreads();
    omp_set_num_threads(availableThreads);

#pragma omp parallel
    {
        int pinOrder = omp_get_thread_num();
        ce->pinThread(pinOrder);
    }
}


