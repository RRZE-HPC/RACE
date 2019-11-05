#include "kernels.h"
#include <mkl.h>
#include "config_eg.h"

inline void SPMV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        double tmp = 0;
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            tmp += mat->val[idx]*x->val[mat->col[idx]];
        }
        b->val[row] = tmp;
    }
}


//b=A*x
void spmv(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    RACE::Interface *ce = mat->ce;

    ENCODE_TO_VOID(mat,b,x);

    int spmvId = ce->registerFunction(&SPMV_KERNEL, voidArg);

    for(int i=0; i<iterations; ++i)
    {
        ce->executeFunction(spmvId);
    }

    DELETE_ARG();
}

inline void SPMTV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        double x_row = x->val[row];
#pragma simd
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            b->val[mat->col[idx]] += mat->val[idx]*x_row;
        }
    }
}


//b=A'*x
void spmtv(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    RACE::Interface *ce = mat->ce;

    ENCODE_TO_VOID(mat,b,x);

    int spmtvId = ce->registerFunction(&SPMTV_KERNEL, voidArg);

    for(int i=0; i<iterations; ++i)
    {
        ce->executeFunction(spmtvId);
    }

    DELETE_ARG();
}

inline void GS_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        x->val[row] = b->val[row];
        double x_row = x->val[row];
        int diag_idx = mat->rowPtr[row];
#pragma simd
        for(int idx=mat->rowPtr[row]+1; idx<mat->rowPtr[row+1]; ++idx)
        {
            x->val[row] -= mat->val[idx]*x->val[mat->col[idx]];
        }
        x->val[row] /= mat->val[diag_idx];
    }
}

//Solve for x : A*x=b
void gs(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    RACE::Interface *ce = mat->ce;

    ENCODE_TO_VOID(mat,b,x);

    int gsId = ce->registerFunction(&GS_KERNEL, voidArg);

    for(int i=0; i<iterations; ++i)
    {
        ce->executeFunction(gsId);
    }

    DELETE_ARG();
}

inline void KACZ_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        double rowNorm = 0.0;
        double scale = 0.0;

        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            double mval = mat->val[idx];
            scale += mval * x->val[mat->col[idx]];
            rowNorm += mval*mval;
        }
        scale /= rowNorm; //omega considered 1

#pragma simd
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            x->val[mat->col[idx]] = x->val[mat->col[idx]] - scale*mat->val[idx];
        }
    }
}

//Solve for x : A*x=b
void kacz(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    RACE::Interface *ce = mat->ce;

    ENCODE_TO_VOID(mat,b,x);

    int kaczId = ce->registerFunction(&KACZ_KERNEL, voidArg);

    for(int i=0; i<iterations; ++i)
    {
        ce->executeFunction(kaczId);
    }

    DELETE_ARG();
}

inline void SYMM_SPMV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        double x_row = x->val[row];
        b->val[row] += mat->val_symm[mat->rowPtr_symm[row]]*x_row;
        double temp = 0;
#pragma simd reduction(+:temp)
        for(int idx=mat->rowPtr_symm[row]+1; idx<mat->rowPtr_symm[row+1]; ++idx)
        {
            double mval = mat->val_symm[idx];
            int colIdx = mat->col_symm[idx];
            temp += mval*x->val[colIdx];
            b->val[colIdx] += mval*x_row;
        }
        b->val[row]+=temp;
   }
}

//A*x=b; A is symmetric
void symm_spmv(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    RACE::Interface *ce = mat->ce;

    ENCODE_TO_VOID(mat,b,x);

    int symm_spmv_Id = ce->registerFunction(&SYMM_SPMV_KERNEL, voidArg);

    for(int i=0; i<iterations; ++i)
    {
        ce->executeFunction(symm_spmv_Id);
    }

    DELETE_ARG();
}

inline void PLAIN_SPMV_KERNEL(int start, int end, int pow, void* args)
{
    DECODE_FROM_VOID(args);

    int parentId = omp_get_thread_num();
//#pragma omp parallel
//    {
        //int ctr = 0;
#pragma omp parallel for schedule(static)
        for(int row=start; row<end; ++row)
        {
            /*
            if(((sched_getcpu())/10) != parentId)
            //if(ctr<100)
            {
                printf("here inside is cpu: %d thread: %d parent: %d\n", sched_getcpu(), omp_get_thread_num(), parentId);
                //++ctr;
            }*/
            double tmp = 0;
            const int offset = (pow-1)*mat->nrows;
#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)
            for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
            {
                tmp += mat->val[idx]*x->val[offset+mat->col[idx]];
            }
            x->val[(pow)*mat->nrows+row] = tmp;
        }
  //  }
}

//plain spmv
void plain_spmv(sparsemat* mat, densemat* x)
{
    ENCODE_TO_VOID(mat, NULL, x);
    PLAIN_SPMV_KERNEL(0, mat->nrows, 1, voidArg);
    DELETE_ARG();
}


inline void MKL_SPMV_KERNEL(int start, int end, int pow, void* args)
{
    DECODE_FROM_VOID(args);
    double *x_ = x->val;
    double *b_ = &(x->val[mat->nrows]);
    int nrows= end-start;

    int nthreads;
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }

    mkl_set_num_threads(nthreads);
    mkl_cspblas_dcsrgemv("N", &nrows, ((double*)mat->val), ((MKL_INT*)mat->rowPtr), ((MKL_INT*) mat->col), ((double*)x_), ((double*)b_));
}

//MKL spmv
void mkl_spmv(sparsemat* mat, densemat* x)
{
    ENCODE_TO_VOID(mat, NULL, x);
    MKL_SPMV_KERNEL(0, mat->nrows, 1, voidArg);
    DELETE_ARG();
}


void matPower(sparsemat *A, int power, densemat *x)
{
    RACE::Interface *ce = A->ce;
    ENCODE_TO_VOID(A, NULL, x);
    int race_power_id = ce->registerFunction(&PLAIN_SPMV_KERNEL, voidArg, power);
//#pragma omp parallel
    {
        ce->executeFunction(race_power_id);
    }
    DELETE_ARG();
}
