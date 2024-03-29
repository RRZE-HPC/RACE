#include "sparsemat.h"
#include "mmio.h"
#include "stdlib.h"
#include <omp.h>
#include <vector>
#include <sys/mman.h>
#include "config_eg.h"
#include <mkl.h>
#ifdef RACE_USE_SPMP
    #include "SpMP/CSR.hpp"
    #include "SpMP/reordering/BFSBipartite.hpp"
#endif
#include "timer.h"
#include "kernels.h"
#include "densemat.h"
#include "multicoloring.h"
#ifdef RACE_HAVE_METIS
    #include "metis.h"
#endif

sparsemat::sparsemat():nrows(0), nnz(0), ce(NULL), val(NULL), rowPtr(NULL), col(NULL), nnz_symm(0), rowPtr_symm(NULL), col_symm(NULL), val_symm(NULL), diagFirst(false), colorType("RACE"), colorBlockSize(64), colorDist(-1), ncolors(-1), colorPtr(NULL), partPtr(NULL), block_size(1), rcmInvPerm(NULL), rcmPerm(NULL), finalPerm(NULL), finalInvPerm(NULL), symm_hint(false)
{
}

//to transfer from a different library the data structure
//need to be called after contructor
void sparsemat::initCover(int nrows_, int nnz_, double *val_, int *rowPtr_, int *col_)
{
    nrows=nrows_;
    nnz=nnz_;
    val=val_;
    rowPtr=rowPtr_;
    col=col_;
}

//performs deep copy of basic data structure
void sparsemat::basicDeepCopy(sparsemat *otherMat)
{
    nrows=otherMat->nrows;
    nnz=otherMat->nnz;

    rowPtr = new int[nrows+1];
    col = new int[nnz];
    val = new double[nnz];

    rowPtr[0] = otherMat->rowPtr[0];
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        rowPtr[row+1] = otherMat->rowPtr[row+1];
        for(int idx=otherMat->rowPtr[row]; idx<otherMat->rowPtr[row+1]; ++idx)
        {
            val[idx] = otherMat->val[idx];
            col[idx] = otherMat->col[idx];
        }
    }
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

    if(finalPerm)
        delete[] finalPerm;

    if(finalInvPerm)
        delete[] finalInvPerm;

/*    if(rowPtr_bcsr)
        delete[] rowPtr_bcsr;

    if(col_bcsr)
        delete[] col_bcsr;

    if(val_bcsr)
        delete[] val_bcsr;*/
}

void sparsemat::printTree()
{
    ce->printZoneTree();
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
    bool pattern_flag = mm_is_pattern(matcode);

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

    //int ncols;
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
        ERROR_PRINT("Matrix not square. Currently only square matrices are supported\n");
        return false;
    }

    //If matrix market file is symmetric; create a general one out of it
    if(symm_flag)
    {
        printf("Creating a general matrix out of a symmetric one\n");
        symm_hint = true;

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
    delete[] perm;

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
    delete[] nnzPerRow;

    numaInit(false);
    //writeFile("beforePerm.mtx");
    return true;
}

bool sparsemat::convertToBCSR(int b_r)
{
    if(rcmPerm)
    {
        permute(rcmPerm, rcmInvPerm);
    }
    bool success_flag = true;
    sparse_matrix_t A_CSR, A_BSR;
    sparse_status_t status;
    //Using MKL to convert to BSR
    status = mkl_sparse_d_create_csr (&A_CSR, SPARSE_INDEX_BASE_ZERO, nrows, nrows, &(rowPtr[0]), &(rowPtr[1]), col, val);
    if(status != SPARSE_STATUS_SUCCESS)
    {
        printf("Error in MKL CRS creation\n");
        success_flag = false;
    }

    status = mkl_sparse_convert_bsr (A_CSR, b_r, SPARSE_LAYOUT_ROW_MAJOR, SPARSE_OPERATION_NON_TRANSPOSE, &A_BSR);
    if(status != SPARSE_STATUS_SUCCESS)
    {
        printf("Error in conversion from CSR to BCSR\n");
        if(status == SPARSE_STATUS_NOT_INITIALIZED)
        {
            printf("Not initialized\n");
        }
        else if(status == SPARSE_STATUS_ALLOC_FAILED)
        {
            printf("Alloc failed\n");
        }
        else if(status == SPARSE_STATUS_INVALID_VALUE)
        {
            printf("invalid value\n");
        }
        else if(status == SPARSE_STATUS_EXECUTION_FAILED)
        {
            printf("execution failed\n");
        }
        else if(status == SPARSE_STATUS_INTERNAL_ERROR)
        {
            printf("Internal error\n");
        }
        else if(status == SPARSE_STATUS_NOT_SUPPORTED)
        {
            printf("Not supported\n");
        }

        success_flag = false;
    }

    status = mkl_sparse_destroy(A_CSR);

    sparse_index_base_t indexing;
    sparse_layout_t block_layout;
    int *rowPtr_start_bcsr_mkl;
    int *rowPtr_end_bcsr_mkl;
    int *col_bcsr_mkl;
    double *val_bcsr_mkl;
    int nrows_bcsr;
    int ncols_bcsr;
    status = mkl_sparse_d_export_bsr (A_BSR, &indexing, &block_layout, &nrows_bcsr, &ncols_bcsr, &block_size, &rowPtr_start_bcsr_mkl, &rowPtr_end_bcsr_mkl, &col_bcsr_mkl, &val_bcsr_mkl);

    if(status != SPARSE_STATUS_SUCCESS)
    {
        printf("Error in exporting to BCSR\n");
        success_flag = false;
    }
    //printf("check rowPtr start = %p, end = %p\n", rowPtr_start_bcsr, rowPtr_end_bcsr);

    if(block_size != b_r)
    {
        printf("Seems like there was an error in CSR to BCSR conversion, requested block size = %d, got %d\n", b_r, block_size);
        success_flag = false;
    }

    if(success_flag)
    {
        nrows = nrows_bcsr;
        int old_nnz = nnz;
        nnz = rowPtr_start_bcsr_mkl[nrows_bcsr];
        printf("Extra nnz for BCSR = %d, i.e. %f\%\n", b_r*b_r*nnz-old_nnz, (b_r*b_r*nnz-old_nnz)*100/(double)old_nnz);

        if(rowPtr != NULL)
        {
            delete[] rowPtr;
        }
        rowPtr = new int[nrows+1];

        for(int i=0; i<nrows+1; ++i)
        {
            rowPtr[i] = rowPtr_start_bcsr_mkl[i];
        }

        if(col != NULL)
        {
            delete[] col;
        }
        if(val != NULL)
        {
            delete[] val;
        }
        col = new int[nnz];
        val = new double[b_r*b_r*nnz];
        for(int i=0; i<nnz; ++i)
        {
            col[i] = col_bcsr_mkl[i];
        }
        for(int i=0; i<b_r*b_r*nnz; ++i)
        {
            val[i] = val_bcsr_mkl[i];
        }
    }

    mkl_sparse_destroy(A_BSR);

    if(success_flag && rcmPerm)
    {
        delete[] rcmPerm;
        delete[] rcmInvPerm;
        rcmPerm = NULL;
        rcmInvPerm = NULL;
        //doRCM();
    }
    return success_flag;
}

bool sparsemat::isAnyDiagZero()
{
    //check whether a new allocation is necessary
    int extra_nnz=0;
    for(int row=0; row<nrows; ++row)
    {
        bool diagHit = false;
        for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
        {
            if(col[idx] == row)
            {
                if(val[idx] != 0)
                {
                    diagHit = true;
                }
            }
        }
        if(!diagHit)
        {
            return true;
        }
    }
    return false;
}

//necessary for GS like kernels
void sparsemat::makeDiagFirst(double missingDiag_value, bool rewriteAllDiag_with_maxRowSum)
{
    double maxRowSum=0.0;
    if(!diagFirst || rewriteAllDiag_with_maxRowSum)
    {
        //check whether a new allocation is necessary
        int extra_nnz=0;
        std::vector<double>* val_with_diag = new std::vector<double>();
        std::vector<int>* col_with_diag = new std::vector<int>();
        std::vector<int>* rowPtr_with_diag = new std::vector<int>(rowPtr, rowPtr+nrows+1);

        for(int row=0; row<nrows; ++row)
        {
            bool diagHit = false;
            double rowSum=0;
            for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
            {
                val_with_diag->push_back(val[idx]);
                col_with_diag->push_back(col[idx]);
                rowSum += val[idx];

                if(col[idx] == row)
                {
                    diagHit = true;
                }
            }
            if(!diagHit)
            {
                val_with_diag->push_back(missingDiag_value);
                col_with_diag->push_back(row);
                ++extra_nnz;
                rowSum += missingDiag_value;
            }
            maxRowSum = std::max(maxRowSum, std::abs(rowSum));
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
                    if(rewriteAllDiag_with_maxRowSum)
                    {
                        newVal[0] = maxRowSum;
                    }
                    else
                    {
                        newVal[0] = val[idx];
                    }
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
        diagFirst = true;
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

            if(row_1_based[idx] < 1)
            {
                printf("row less than 1: value=%d at row=%d, idx=%d\n", i+1, i, idx);
            }
            if(col_1_based[idx] < 1)
            {
                printf("col less than 1: value=%d at row=%d,idx=%d\n", col[idx]+1, i, idx);
            }

        }
    }

    mm_write_mtx_crd(filename, nrows, nrows, nnz, row_1_based, col_1_based, val, "MCRG");

    delete[] row_1_based;
    delete[] col_1_based;

    return true;
}

bool sparsemat::computeSymmData()
{
    //this is assumed by SymmSpMV kernel and GS
    makeDiagFirst();

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


void sparsemat::splitMatrixToLU(sparsemat **L_ptr, sparsemat **U_ptr)
{

    int* L_rowPtr = new int[nrows+1];
    int* U_rowPtr = new int[nrows+1];

    L_rowPtr[0] = 0;
    U_rowPtr[0] = 0;

    //NUMA init
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        L_rowPtr[row+1] = 0;
        U_rowPtr[row+1] = 0;
    }

    int L_nnz = 0;
    int U_nnz = 0;
    for(int row=0; row<nrows; ++row)
    {
        for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
        {
            if(col[idx] > row)
            {
                ++U_nnz;
            }
            else
            {
                ++L_nnz;
            }
        }
        L_rowPtr[row+1] = L_nnz;
        U_rowPtr[row+1] = U_nnz;
    }

    double* L_val = new double[L_nnz];
    int* L_col = new int[L_nnz];
    double* U_val = new double[U_nnz];
    int* U_col = new int[U_nnz];

    //with NUMA init
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        int L_ctr = L_rowPtr[row];
        int U_ctr = U_rowPtr[row];
        for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
        {
            if(col[idx]>row)
            {
                U_col[U_ctr] = col[idx];
                U_val[U_ctr] = val[idx];
                ++U_ctr;
            }
            else
            {
                L_col[L_ctr] = col[idx];
                L_val[L_ctr] = val[idx];
                ++L_ctr;
            }
        }
    }

    (*L_ptr)->initCover(nrows, L_nnz, L_val, L_rowPtr, L_col);
    (*U_ptr)->initCover(nrows, U_nnz, U_val, U_rowPtr, U_col);

}

bool sparsemat::isSymmetric()
{
#ifdef RACE_USE_SPMP
    SpMP::CSR *csr = NULL;
    csr = new SpMP::CSR(nrows, nrows, rowPtr, col, val);
    return csr->isSymmetric(true, true);
#else
    printf("Please link with SpMP library to check for symmetry.\n");
    return -1;
#endif
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
 //   rcmPerm = new int[nrows];
 //    rcmInvPerm = new int[nrows];
    if(csr->isSymmetric(false,false))
    {
        rcmPerm = new int[nrows];
        rcmInvPerm = new int[nrows];
        csr->getRCMPermutation(rcmInvPerm, rcmPerm);
    }
    else
    {
        printf("Matrix not symmetric RCM cannot be done\n");
    }
    omp_set_num_threads(orig_threads);
    delete csr;
#else
    printf("Please link with SpMP library to enable RCM permutation\n");
#endif
}

void sparsemat::doRCMPermute()
{
    doRCM();
    permute(rcmPerm, rcmInvPerm);
    if(finalPerm)
    {
        delete [] finalPerm;
        delete [] finalInvPerm;
    }
    finalPerm = rcmPerm;
    finalInvPerm = rcmInvPerm;
    rcmPerm = NULL;
    rcmInvPerm = NULL;
}

int sparsemat::prepareForPower(int highestPower, double cacheSize, int nthreads, int smt, PinMethod pinMethod, int globalStartRow, int globalEndRow, std::string mtxType)
{
    //permute(rcmInvPerm, rcmPerm);
    //rcmPerm = NULL;
    //rcmInvPerm = NULL;
    INIT_TIMER(pre_process_kernel);
    START_TIMER(pre_process_kernel);
    ce = new Interface(nrows, nthreads, RACE::POWER, rowPtr, col, symm_hint, smt, pinMethod, rcmPerm, rcmInvPerm);
    //ce->RACEColor(highestPower, cacheSize);
    ce->RACEColor(highestPower, cacheSize*1024*1024, 2, mtxType);
    //ce->RACEColor(highestPower, cacheSize*1024*1024, 2, mtxType, 3);
    if ((globalStartRow != -1) && (globalEndRow != -1))
        ce->passGlobalRows(globalStartRow, globalEndRow);
    STOP_TIMER(pre_process_kernel);
    printf("Pre-processing time: cache size = %f, power = %d, RACE pre-processing time = %fs\n", cacheSize, highestPower, GET_TIMER(pre_process_kernel));

    int *perm, *invPerm, permLen;
    ce->getPerm(&perm, &permLen);
    ce->getInvPerm(&invPerm, &permLen);
    permute(perm, invPerm, true);

    if(finalPerm)
    {
        delete [] finalPerm;
        delete [] finalInvPerm;
    }

    finalPerm = perm;
    finalInvPerm = invPerm;

    checkNumVecAccesses(highestPower);
    //delete [] invPerm;
    //delete [] perm;
    //no idea why need it second time w/o perm. 
    //NUMA init work nicely only if this is done; (only for pwtk, others perf
    //degradation))
    //numaInit(true);
    //writeFile("after_RCM.mtx");
    return 1;
}

struct metis_block_arg
{
    int blockSize;
    int* init_perm;
    int* init_invPerm;
    int* out_invPerm;
    sparsemat* mat;
};

inline void METIS_BLOCKING_KERNEL(int start, int end, void *args)
{
    metis_block_arg* arg_decode = (metis_block_arg*) args;
    if(start<end)
    {
        printf("start=%d, end=%d\n", start, end);
        arg_decode->mat->doMETIS(arg_decode->blockSize, start, end, arg_decode->init_perm, arg_decode->init_invPerm, arg_decode->out_invPerm);
    }
}

int sparsemat::colorAndPermute(dist distance, std::string colorType_, int nthreads, int smt, PinMethod pinMethod)
{
    INIT_TIMER(pre_process_kernel);
    START_TIMER(pre_process_kernel);

    //writeFile("before_perm.mtx");
    colorType = colorType_;
    if((colorType == "RACE") || (colorType == "RACE_and_BLOCK"))
    {
        ce = new Interface(nrows, nthreads, distance, rowPtr, col, symm_hint, smt, pinMethod, rcmPerm, rcmInvPerm);
        RACE_error ret = ce->RACEColor();

        if(ret != RACE_SUCCESS)
        {
            printf("Pinning failure\n");
            return 0;
        }

        int *perm, *invPerm, permLen;

        ce->getPerm(&perm, &permLen);
        ce->getInvPerm(&invPerm, &permLen);

        if(finalPerm)
        {
            delete [] finalPerm;
            delete [] finalInvPerm;
        }

        finalPerm = perm;
        finalInvPerm = invPerm;

        if(colorType == "RACE_and_BLOCK")
        {
            int* metis_perm=new int[nrows];
            int* metis_invPerm=new int[nrows];
            int blockSize=128;

            metis_block_arg* arg_encode = new metis_block_arg;
            arg_encode->blockSize=blockSize;
            arg_encode->init_perm=perm;
            arg_encode->init_invPerm=invPerm;
            arg_encode->out_invPerm=metis_invPerm;
            arg_encode->mat=this;
            void *voidArg = (void*) arg_encode;

            int metis_block_kernel_id = ce->registerFunction(&METIS_BLOCKING_KERNEL, voidArg);
            ce->executeFunction(metis_block_kernel_id);
            for (int i=0; i<nrows; i++)
            {
                metis_perm[metis_invPerm[i]] = i;
            }

            if(finalPerm)
            {
                delete [] finalPerm;
                delete [] finalInvPerm;
            }

            finalPerm = metis_perm;
            finalInvPerm = metis_invPerm;
        }

        //pin omp threads as in RACE for proper NUMA
        //pinOMP(nthreads);

        double parallel_efficiency = ce->getEfficiency();
        printf("Parallel efficiency of RACE = %f\n", parallel_efficiency);


    }
    else
    {
        int dist_int = 1;
        if(distance == RACE::ONE)
        {
            dist_int = 1;
        }
        else
        {
            dist_int = 2;
        }
        multicoloring mc(this, dist_int, rcmPerm, rcmInvPerm);

        if(colorType == "MC")
        {
            mc.doMC();
        }
        else
        {
            mc.doABMC();
            partPtr = mc.partPtr_out;
        }

        ncolors = mc.ncolors_out;
        colorPtr = mc.colorPtr_out;
/*
        for(int color=0; color<ncolors+1; ++color)
        {
            printf("check colorPtr[%d] = %d\n", color, colorPtr[color]);
        }

        for(int row=0; row<nrows; ++row)
        {
            printf("check perm[%d] = %d\n", row, mc.perm_out[row]);
        }
*/

        if(finalPerm)
        {
            delete[] finalPerm;
            delete[] finalInvPerm;
        }

        finalPerm = mc.perm_out;
        finalInvPerm = mc.invPerm_out;
    }

    STOP_TIMER(pre_process_kernel);
    printf("RACE pre-processing time = %fs\n", GET_TIMER(pre_process_kernel));

    permute(finalPerm, finalInvPerm);

    //writeFile("after_perm.mtx");
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

densemat* sparsemat::permute_densemat(densemat *vec)
{
    densemat* newVec = new densemat(nrows, vec->ncols);
    newVec->setVal(0);

    if(nrows != vec->nrows)
    {
        ERROR_PRINT("Permutation of densemat not possible, dimension of matrix and vector do not match");
    }
    for(int row=0; row<nrows; ++row)
    {
        int perm_row = row;
        if(finalPerm != NULL)
        {
            perm_row = finalPerm[row];
        }
        for(int col=0; col<vec->ncols; ++col)
        {
            newVec->val[col*nrows+row] = vec->val[col*nrows+perm_row];
        }
    }
    return newVec;
}

//symmetrically permute
void sparsemat::permute(int *_perm_, int*  _invPerm_, bool RACEalloc)
{
    double* newVal = (double*)malloc(sizeof(double)*block_size*block_size*nnz);
        //new double[block_size*block_size*nnz];
    int* newCol = (int*)malloc(sizeof(int)*nnz);
        //new int[nnz];
    int* newRowPtr = (int*)malloc(sizeof(int)*(nrows+1));
        //new int[nrows+1];

/*
    double *newVal = (double*) malloc(sizeof(double)*nnz);
    int *newCol = (int*) malloc(sizeof(int)*nnz);
    int *newRowPtr = (int*) malloc(sizeof(int)*(nrows+1));
*/

    newRowPtr[0] = 0;

    printf("Initing rowPtr\n");
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

    if(_perm_ != NULL)
    {
        //first find newRowPtr; therefore we can do proper NUMA init
        int _perm_Idx=0;
        printf("nrows = %d\n", nrows);
        for(int row=0; row<nrows; ++row)
        {
            //row _perm_utation
            int _perm_Row = _perm_[row];
            for(int idx=rowPtr[_perm_Row]; idx<rowPtr[_perm_Row+1]; ++idx)
            {
                ++_perm_Idx;
            }
            newRowPtr[row+1] = _perm_Idx;
        }
    }
    else
    {
        for(int row=0; row<nrows+1; ++row)
        {
            newRowPtr[row] = rowPtr[row];
        }
    }


    printf("Initing mtxVec\n");
    if(RACEalloc)
    {
        ce->numaInitMtxVec(newRowPtr, newCol, newVal, NULL);
    }

    printf("Finished inting\n");
    if(_perm_ != NULL)
    {
        //with NUMA init
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            //row _perm_utation
            int _perm_Row = _perm_[row];
            for(int _perm_Idx=newRowPtr[row],idx=rowPtr[_perm_Row]; _perm_Idx<newRowPtr[row+1]; ++idx,++_perm_Idx)
            {
                //_perm_ute column-wise also
                //newVal[_perm_Idx] = val[idx];
                newCol[_perm_Idx] = _invPerm_[col[idx]];
                for(int b=0; b<block_size*block_size; ++b)
                {
                    newVal[_perm_Idx+b] = val[idx+b];
                }
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
                //newVal[idx] = val[idx];
                newCol[idx] = col[idx];
                for(int b=0; b<block_size*block_size; ++b)
                {
                    newVal[idx+b] = val[idx+b];
                }
            }
        }
    }


    //free old _perm_utations
    delete[] val;
    delete[] rowPtr;
    delete[] col;

    val = newVal;
    rowPtr = newRowPtr;
    col = newCol;
}

//make sure all permutation are actuakly done before entring this routine,
//because no initPerm is supported now
//also outPerm and outInvPerm have to be allocated before entering this
//routine
bool sparsemat::doMETIS(int blocksize, int start_row, int end_row, int *initPerm, int *initInvPerm, int *outInvPerm)
{
#if (defined RACE_HAVE_METIS)
    bool success_flag=true;
    int* rptlocal = (int*) malloc(sizeof(int)*(nrows+1));
    int* collocal = (int*) malloc(sizeof(int)*nnz);
    rptlocal[0] = 0;
    int local_nrows = end_row-start_row;
    int local_nnz = 0;
    int local_row=0;

    for(int row=start_row; row<end_row; ++row, ++local_row)
    {
        int permRow = row;
        if(initPerm)
        {
            permRow = initPerm[row];
        }
        for(int idx=rowPtr[permRow]; idx<rowPtr[permRow+1]; ++idx)
        {
            int permCol = col[idx];
            if(initInvPerm)
            {
                permCol = initInvPerm[permCol];
            }
            //only local entries added to col
            if((permCol>=start_row) && (permCol<end_row))
            {
                collocal[local_nnz] = permCol-start_row; //-start_row to convert to local index
                ++local_nnz;
            }
        }
        rptlocal[local_row+1] = local_nnz;
    }

    //partition using METIS
    int ncon = 1;
    int nparts = (int)(local_nrows/(double)blocksize);
    int objval;
    int *part = (int*) malloc(sizeof(int)*local_nrows);

    printf("partitioning graph to %d parts\n", nparts);
    int metis_ret = METIS_PartGraphKway(&local_nrows, &ncon, rptlocal, collocal, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);
    if(metis_ret == METIS_OK)
    {
        printf("successfully partitioned graph to nparts=%d\n", nparts);
    }
    else
    {
        success_flag=false;
        printf("Error in METIS partitioning\n");
    }

    std::vector<std::vector<int>> partRow(nparts);
    for (int i=0;i<local_nrows;i++) {
        partRow[part[i]].push_back(i);
    }
/*    for(int i=start_row; i<end_row; ++i) {
        outPerm[i]=-1;
        outInvPerm[i]=-1;
    }*/
    int ctr=0;
    for (int partIdx=0; partIdx<nparts; partIdx++)
    {
        int partSize = (int)partRow[partIdx].size(); //partPtr[partIdx+1]-partPtr[partIdx];
        for(int rowIdx=0; rowIdx<partSize; ++rowIdx)
        {
            //find rows in parts
            int local_currRow = partRow[partIdx][rowIdx];
            int currRow = local_currRow + start_row;
            int permCurrRow = currRow;
            if(initPerm)
            {
                permCurrRow = initPerm[currRow];
            }

            outInvPerm[permCurrRow] = start_row+ctr;
            ++ctr;
        }
    }
/*
    for (int i=start_row; i<end_row; i++)
    {
        outPerm[outInvPerm[i]] = i;
    }
*/
    return success_flag;
#else
    return false;
#endif
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


NUMAmat::NUMAmat(sparsemat *mat_, bool manual, std::vector<int> splitRows_):mat(mat_)
{
    if(manual)
    {
        splitRows = splitRows_;
    }
    else
    {
        splitRows = getRACEPowerSplit();
    }
    NUMAdomains  = splitRows.size()-1;

    nrows = new int[NUMAdomains];
    nnz = new int[NUMAdomains];
    rowPtr = new int*[NUMAdomains];
    col = new int*[NUMAdomains];
    val = new double*[NUMAdomains];
    for(int domain = 0; domain<NUMAdomains; ++domain)
    {
        int startRow = splitRows[domain];
        int endRow = splitRows[domain+1];
        int currNrows = endRow - startRow;
        rowPtr[domain] = new int[currNrows+1];
    }

    if(!manual)
    {
    //    mat->ce->numaInitRowPtr(rowPtr);
    }

#pragma omp parallel
    {
        int totalThreads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        int threadPerNode = totalThreads/NUMAdomains;
        int domain = tid/threadPerNode;
        int localTid = tid%threadPerNode;

        if(localTid == 0)
        {
            int startRow = splitRows[domain];
            int endRow = splitRows[domain+1];
            int currNrows = endRow - startRow;
            //BCSR not yet for NUMAmat
            nrows[domain] = currNrows;
            //rowPtr[domain] = new int[currNrows+1];
            rowPtr[domain][0] = 0;
            int cur_nnz = 0;
            for(int row=0; row<currNrows; ++row)
            {
                //rowPtr[domain][row+1] = mat->rowPtr[row+1+startRow];

                for(int idx=mat->rowPtr[row+startRow]; idx<mat->rowPtr[row+1+startRow]; ++idx)
                {
                    ++cur_nnz;
                }
                rowPtr[domain][row+1] = cur_nnz;
            }
            nnz[domain] = cur_nnz;
        }
    }
    for(int domain = 0; domain<NUMAdomains; ++domain)
    {
        int cur_nnz = nnz[domain];
        col[domain] = new int[cur_nnz];
        val[domain] = new double[cur_nnz];
    }

    if(!manual)
    {
    //    mat->ce->numaInitMtxVec(rowPtr, col, val, 1);
    }

#pragma omp parallel
    {
        int totalThreads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        int threadPerNode = totalThreads/NUMAdomains;
        int domain = tid/threadPerNode;
        int localTid = tid%threadPerNode;

        if(localTid == 0)
        {
            int startRow = splitRows[domain];
            int endRow = splitRows[domain+1];
            int currNrows = endRow - startRow;

            for(int row=0; row<currNrows; ++row)
            {
                for(int idx=mat->rowPtr[row+startRow], local_idx=rowPtr[domain][row]; idx<mat->rowPtr[row+1+startRow]; ++idx,++local_idx)
                {
                    col[domain][local_idx] = mat->col[idx];
                    val[domain][local_idx] = mat->val[idx];
                }
            }
        }
    }
    /*
    printf("Mat = \n");
    for(int row=0; row<mat->nrows; ++row)
    {
        printf("row = %d \t", row);
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            printf("(%d,%f) ", mat->col[idx], mat->val[idx]);
        }
        printf("\n");
    }

    printf("Sub mat = \n");
    for(int domain=0; domain<NUMAdomains; ++domain)
    {
        printf("domain = %d\n", domain);
        for(int row=0; row<nrows[domain]; ++row)
        {
            printf("row = %d \t", row);
            for(int idx=rowPtr[domain][row]; idx<rowPtr[domain][row+1]; ++idx)
            {
                printf("(%d,%f) ", col[domain][idx], val[domain][idx]);
            }
            printf("\n");
        }
    }
*/
}

std::vector<int> NUMAmat::getRACEPowerSplit()
{
    int *split, splitLen;
    mat->ce->getNumaSplitting(&split, &splitLen);
    std::vector<int> split_vec;
    for(int i=0; i<splitLen; ++i)
    {
        split_vec.push_back(split[i]);
    }

    return split_vec;
}

NUMAmat::~NUMAmat()
{
    for(int domain = 0; domain<NUMAdomains; ++domain)
    {
        if(rowPtr[domain])
        {
            delete[] rowPtr[domain];
        }

        if(col[domain])
        {
            delete[] col[domain];
        }

        if(val[domain])
        {
            delete[] val[domain];
        }
    }
    delete[] rowPtr;
    delete[] col;
    delete[] val;
}

inline void MAT_NUM_VEC_ACCESSES(int start, int end, int pow, int numa_domain, void* args)
{
    DECODE_FROM_VOID(args);

    int nrows=mat->nrows;
    for(int row=start; row<end; ++row)
    {
        x->val[row]++;
        if(x->val[row] != pow)
        {
            if(x->val[row] > pow)
            {
                ERROR_PRINT("Oh oh we have duplicate computations, error at pow=%d, for row=%d. Value I got at x is %f, expected %d. Level start =%d, Level end=%d", pow, row, x->val[row], pow, start, end);
            }
            else
            {
                ERROR_PRINT("Oh oh have some missing computations, error at pow=%d, for row=%d. Value I got at x is %f, expected %d. Level start =%d, Level end=%d", pow, row, x->val[row], pow, start, end);
            }
        }
    }
}

inline void MAT_NUM_VEC_ACCESSES_w_subPower(int start, int end, int pow, int subPow, int numa_domain, void* args)
{
    DECODE_FROM_VOID(args);

    int nrows=mat->nrows;
    for(int row=start; row<end; ++row)
    {
        x->val[row]++;
        if(x->val[row] != (pow-1)*3+subPow)
        {
            ERROR_PRINT("Oh oh we have duplicate computations, error at pow=%d, subPow=%d, for row=%d. Value I got at x is %f, expected %d. Level start =%d, Level end=%d", pow, subPow, row, x->val[row], (pow-1)*3+subPow, start, end);
        }
    }
}

void sparsemat::checkNumVecAccesses(int power)
{
    densemat* x = new densemat(this->nrows);
    ENCODE_TO_VOID(this, NULL, x);
    int race_power_id = ce->registerFunction(&MAT_NUM_VEC_ACCESSES, voidArg, power);
    //int race_power_id = ce->registerFunction(&MAT_NUM_VEC_ACCESSES_w_subPower, voidArg, power, 3);
    {
        ce->executeFunction(race_power_id);
    }
    DELETE_ARG();
    delete x;
}


