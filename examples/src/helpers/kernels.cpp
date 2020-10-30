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

void plain_spmv(densemat* b, sparsemat* mat, densemat* x, int iterations)
{

    for(int i=0; i<iterations; ++i)
    {
#pragma omp parallel for schedule(static)
        for(int row=0; row<mat->nrows; ++row)
        {
            double tmp = 0;
            for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
            {
                tmp += mat->val[idx]*x->val[mat->col[idx]];
            }
            b->val[row] = tmp;
        }
    }
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

        scale = -b->val[row];
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


inline void PLAIN_SPMV_NUMA_KERNEL(int start, int end, int pow, void* args)
{
    DECODE_FROM_VOID_NUMA(args);

    int parentId = omp_get_thread_num();
//#pragma omp parallel
//    {
        //int ctr = 0;
#pragma omp parallel for schedule(static)
        for(int row=start; row<end; ++row)
        {
            int nrows = mat->mat->nrows;
            int totalThreads = omp_get_num_threads();
            int threadPerNode = totalThreads/mat->NUMAdomains;
            int tid = omp_get_thread_num();
            int numa_domain = tid/threadPerNode;
            double tmp = 0;
            const int col_offset = (pow-1)*nrows;
            const int row_offset = mat->splitRows[numa_domain];
#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)
            for(int idx=mat->rowPtr[numa_domain][row]; idx<mat->rowPtr[numa_domain][row+1]; ++idx)
            {
                tmp += mat->val[numa_domain][idx]*x->val[col_offset+mat->col[numa_domain][idx]];
            }
            x->val[(pow)*nrows+row+row_offset] = tmp;
        }
  //  }
}

void plain_spmv_numa(NUMAmat *mat, densemat *x)
{
    ENCODE_TO_VOID_NUMA(mat, NULL, x);
    PLAIN_SPMV_NUMA_KERNEL(0, mat->mat->nrows, 1, voidArg);
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

inline void MAT_SPMV_KERNEL(int start, int end, int pow, int numa_domain, void* args)
{
    DECODE_FROM_VOID(args);

    int parentId = omp_get_thread_num();
//#pragma omp parallel
//    {
        //int ctr = 0;
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

void matPower(sparsemat *A, int power, densemat *x)
{
    RACE::Interface *ce = A->ce;
    ENCODE_TO_VOID(A, NULL, x);
    int race_power_id = ce->registerFunction(&MAT_SPMV_KERNEL, voidArg, power);
    {
        ce->executeFunction(race_power_id);
    }
    DELETE_ARG();
}

inline void MAT_SPMV_NUMA_KERNEL(int start, int end, int pow, int numa_domain, void* args)
{
    DECODE_FROM_VOID_NUMA(args);

    //int parentId = omp_get_thread_num();
//#pragma omp parallel
//    {
        //int ctr = 0;
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
            int nrows = mat->mat->nrows;
            const int col_offset = (pow-1)*nrows;
            const int row_offset = mat->splitRows[numa_domain];
#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)
            for(int idx=mat->rowPtr[numa_domain][row]; idx<mat->rowPtr[numa_domain][row+1]; ++idx)
            {
                tmp += mat->val[numa_domain][idx]*x->val[col_offset+mat->col[numa_domain][idx]];
            }
            x->val[(pow)*nrows+row+row_offset] = tmp;
        }
  //  }
}

void matPowerNuma(NUMAmat *A, int power, densemat *x)
{
    RACE::Interface *ce = A->mat->ce;
    ENCODE_TO_VOID_NUMA(A, NULL, x);
    int race_power_id = ce->registerFunction(&MAT_SPMV_NUMA_KERNEL, voidArg, power, true);
    {
        ce->executeFunction(race_power_id);
    }
    DELETE_ARG();
}



inline void MAT_SPMV_KERNEL_BCSR_2x2(int start, int end, int pow, int numa_domain, void* args)
{
    DECODE_FROM_VOID(args);

    int parentId = omp_get_thread_num();
    //int ctr = 0;
    for(int row=start; row<end; ++row)
    {
        /*
           if(((sched_getcpu())/10) != parentId)
        //if(ctr<100)
        {
        printf("here inside is cpu: %d thread: %d parent: %d\n", sched_getcpu(), omp_get_thread_num(), parentId);
        //++ctr;
        }*/
        double tmp[2];
        tmp[0] = 0;
        tmp[1] = 0;
        const int offset = (pow-1)*mat->nrows*2;
//#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH)
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
#pragma unroll_and_jam(4)
            for(int b_x=0; b_x<2; ++b_x)
            {
                for(int b_y=0; b_y<2; ++b_y)
                {
                    tmp[b_x] += mat->val[4*idx+b_x*2+b_y]*x->val[offset+2*mat->col[idx]+b_y];
                }
            }
        }
        x->val[(pow)*mat->nrows*2+2*row] = tmp[0];
        x->val[(pow)*mat->nrows*2+2*row+1] = tmp[1];
    }
}

void matPowerBCSR(sparsemat *A, int power, densemat *x)
{
    RACE::Interface *ce = A->ce;
    ENCODE_TO_VOID(A, NULL, x);
    int race_power_id = ce->registerFunction(&MAT_SPMV_KERNEL_BCSR_2x2, voidArg, power);
    {
        ce->executeFunction(race_power_id);
    }
    DELETE_ARG();
}

