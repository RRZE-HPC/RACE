#include "kernels.h"
#include <mkl.h>
#include "config_eg.h"

#define SPMV_KERNEL_BODY()\
{\
    double tmp = 0;\
_Pragma("nounroll")\
_Pragma("simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)")\
    for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)\
    {\
        tmp += mat->val[idx]*x->val[mat->col[idx]];\
    }\
    b->val[row] = tmp;\
}\

inline void SPMV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        SPMV_KERNEL_BODY();
    }
}


//b=A*x
void spmv(densemat* b, sparsemat* mat, densemat* x)
{
    if(mat->colorType == "RACE")
    {
        RACE::Interface *ce = mat->ce;

        ENCODE_TO_VOID(mat,b,x);

        int spmvId = ce->registerFunction(&SPMV_KERNEL, voidArg);

        ce->executeFunction(spmvId);

        DELETE_ARG();
    }
    else if(mat->colorType == "ABMC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int part=mat->colorPtr[color]; part<mat->colorPtr[color+1]; ++part)
            {
                for(int row=mat->partPtr[part]; row<mat->partPtr[part+1]; ++row)
                {
                    SPMV_KERNEL_BODY();
                }
            }
        }
    }
    else if(mat->colorType == "MC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int row=mat->colorPtr[color]; row<mat->colorPtr[color+1]; ++row)
            {
                SPMV_KERNEL_BODY();
            }
        }
    }
    else
    {
        printf("Didn't recognize the coloring type. Available options: RACE, ABMC, MC\n");
    }
}

/*
void plain_spmv(densemat* b, sparsemat* mat, densemat* x)
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
*/

void plain_spmv(densemat* b, sparsemat* mat, densemat* x)
{

#pragma omp parallel for schedule(static)
        for(int row=0; row<mat->nrows; ++row)
        {
            double tmp = 0;
#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)
            for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
            {
                tmp += mat->val[idx]*x->val[mat->col[idx]];
            }
            b->val[row] = tmp;
        }
}

#define SPMTV_KERNEL_BODY()\
{\
    double x_row = x->val[row];\
    _Pragma("simd vectorlength(VECTOR_LENGTH)")\
    for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)\
    {\
        b->val[mat->col[idx]] += mat->val[idx]*x_row;\
    }\
}\

inline void SPMTV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        SPMTV_KERNEL_BODY();
    }
}


//b=A'*x
void spmtv(densemat* b, sparsemat* mat, densemat* x)
{
    if(mat->colorType == "RACE")
    {
        RACE::Interface *ce = mat->ce;

        ENCODE_TO_VOID(mat,b,x);

        int spmtvId = ce->registerFunction(&SPMTV_KERNEL, voidArg);

        ce->executeFunction(spmtvId);

        DELETE_ARG();
    }
    else if(mat->colorType == "ABMC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int part=mat->colorPtr[color]; part<mat->colorPtr[color+1]; ++part)
            {
                for(int row=mat->partPtr[part]; row<mat->partPtr[part+1]; ++row)
                {
                    SPMTV_KERNEL_BODY();
                }
            }
        }
    }
    else if(mat->colorType == "MC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int row=mat->colorPtr[color]; row<mat->colorPtr[color+1]; ++row)
            {
                SPMTV_KERNEL_BODY();
            }
        }
    }
    else
    {
        printf("Didn't recognize the coloring type. Available options: RACE, ABMC, MC\n");
    }
}

#define GS_KERNEL_BODY()\
{\
    double x_row = x->val[row];\
    int diag_idx = mat->rowPtr[row];\
    double tmp=b->val[row];\
_Pragma("simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)")\
    for(int idx=mat->rowPtr[row]+1; idx<mat->rowPtr[row+1]; ++idx)\
    {\
        tmp += (-mat->val[idx])*x->val[mat->col[idx]];\
    }\
    x->val[row] = tmp/mat->val[diag_idx];\
}\

inline void GS_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        GS_KERNEL_BODY();
    }
}

//Solve for x : A*x=b
void gs_serial(densemat* b, sparsemat* mat, densemat* x)
{
    for(int row=0; row<mat->nrows; ++row)
    {
        GS_KERNEL_BODY();
    }
}

//Solve for x : A*x=b
void gs(densemat* b, sparsemat* mat, densemat* x)
{
    if(mat->colorType == "RACE")
    {
        RACE::Interface *ce = mat->ce;

        ENCODE_TO_VOID(mat,b,x);

        int gsId = ce->registerFunction(&GS_KERNEL, voidArg);

        ce->executeFunction(gsId);

        DELETE_ARG();
    }
    else if(mat->colorType == "ABMC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int part=mat->colorPtr[color]; part<mat->colorPtr[color+1]; ++part)
            {
                for(int row=mat->partPtr[part]; row<mat->partPtr[part+1]; ++row)
                {
                    GS_KERNEL_BODY();
                }
            }
        }
    }
    else if(mat->colorType == "MC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int row=mat->colorPtr[color]; row<mat->colorPtr[color+1]; ++row)
            {
                GS_KERNEL_BODY();
            }
        }
    }
    else
    {
        printf("Didn't recognize the coloring type. Available options: RACE, ABMC, MC\n");
    }
}

#define KACZ_KERNEL_BODY()\
{\
    double rowNorm = 0.0;\
    double scale = 0.0;\
    scale = -b->val[row];\
    for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)\
    {\
        double mval = mat->val[idx];\
        scale += mval * x->val[mat->col[idx]];\
        rowNorm += mval*mval;\
    }\
    scale /= rowNorm; /*omega considered 1*/\
_Pragma("simd vectorlength(VECTOR_LENGTH)")\
    for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)\
    {\
        x->val[mat->col[idx]] = x->val[mat->col[idx]] - scale*mat->val[idx];\
    }\
}\

inline void KACZ_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        KACZ_KERNEL_BODY();
    }
}

//Solve for x : A*x=b
void kacz_serial(densemat* b, sparsemat* mat, densemat* x)
{
    for(int row=0; row<mat->nrows; ++row)
    {
        KACZ_KERNEL_BODY();
    }
}

//Solve for x : A*x=b
void kacz(densemat* b, sparsemat* mat, densemat* x)
{
    if(mat->colorType == "RACE")
    {
        RACE::Interface *ce = mat->ce;

        ENCODE_TO_VOID(mat,b,x);

        int kaczId = ce->registerFunction(&KACZ_KERNEL, voidArg);

        ce->executeFunction(kaczId);

        DELETE_ARG();
    }
    else if(mat->colorType == "ABMC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int part=mat->colorPtr[color]; part<mat->colorPtr[color+1]; ++part)
            {
                for(int row=mat->partPtr[part]; row<mat->partPtr[part+1]; ++row)
                {
                    KACZ_KERNEL_BODY();
                }
            }
        }
    }
    else if(mat->colorType == "MC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int row=mat->colorPtr[color]; row<mat->colorPtr[color+1]; ++row)
            {
                KACZ_KERNEL_BODY();
            }
        }
    }
    else
    {
        printf("Didn't recognize the coloring type. Available options: RACE, ABMC, MC\n");
    }
}

#define SYMM_SPMV_KERNEL_BODY()\
{\
    double x_row = x->val[row];\
    b->val[row] += mat->val_symm[mat->rowPtr_symm[row]]*x_row;\
    double temp = 0;\
_Pragma("simd reduction(+:temp)") /*Not doiing here for vector length because nnzr can be very small, lat compiler decide*/\
    for(int idx=mat->rowPtr_symm[row]+1; idx<mat->rowPtr_symm[row+1]; ++idx)\
    {\
        double mval = mat->val_symm[idx];\
        int colIdx = mat->col_symm[idx];\
        temp += mval*x->val[colIdx];\
        b->val[colIdx] += mval*x_row;\
    }\
    b->val[row]+=temp;\
}\

inline void SYMM_SPMV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        SYMM_SPMV_KERNEL_BODY();
    }
}

//A*x=b; A is symmetric
void symm_spmv(densemat* b, sparsemat* mat, densemat* x)
{

    if(mat->colorType == "RACE")
    {
        RACE::Interface *ce = mat->ce;

        ENCODE_TO_VOID(mat,b,x);

        int symm_spmv_Id = ce->registerFunction(&SYMM_SPMV_KERNEL, voidArg);

        ce->executeFunction(symm_spmv_Id);

        DELETE_ARG();
    }
    else if(mat->colorType == "ABMC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int part=mat->colorPtr[color]; part<mat->colorPtr[color+1]; ++part)
            {
                for(int row=mat->partPtr[part]; row<mat->partPtr[part+1]; ++row)
                {
                    SYMM_SPMV_KERNEL_BODY();
                }
            }
        }
    }
    else if(mat->colorType == "MC")
    {
        for(int color=0; color<mat->ncolors; ++color)
        {
#pragma omp parallel for schedule(static)
            for(int row=mat->colorPtr[color]; row<mat->colorPtr[color+1]; ++row)
            {
                SYMM_SPMV_KERNEL_BODY();
            }
        }
    }
    else
    {
        printf("Didn't recognize the coloring type. Available options: RACE, ABMC, MC\n");
    }
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
void MAT_SPMV_KERNEL(int start, int end, int pow, int numa_domain, void* args);

//plain spmv
void plain_spmv(sparsemat* mat, densemat* x)
{
    ENCODE_TO_VOID(mat, NULL, x);
    PLAIN_SPMV_KERNEL(0, mat->nrows, 1, voidArg);
    //MAT_SPMV_KERNEL(0, mat->nrows, 1, 0, voidArg);
    DELETE_ARG();
}

inline void PLAIN_SPMV_KERNEL_only_highest(int start, int end, int pow, void* args)
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
            const int offset = ((pow-1)%2)*mat->nrows;
#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)
            for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
            {
                tmp += mat->val[idx]*x->val[offset+mat->col[idx]];
            }
            x->val[(pow%2)*mat->nrows+row] = tmp;
        }
  //  }
}

//plain spmv
void plain_spmv_only_highest(sparsemat* mat, densemat* x, int pow)
{
    ENCODE_TO_VOID(mat, NULL, x);
    PLAIN_SPMV_KERNEL_only_highest(0, mat->nrows, pow, voidArg);
    //MAT_SPMV_KERNEL(0, mat->nrows, 1, 0, voidArg);
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
        for(int row_global=start; row_global<end; ++row_global)
        {
            int nrows = mat->mat->nrows;
            int totalThreads = omp_get_num_threads();
            int threadPerNode = totalThreads/mat->NUMAdomains;
            int tid = omp_get_thread_num();
            int numa_domain = tid/threadPerNode;

            const int row_offset = mat->splitRows[numa_domain];
            int row = row_global - row_offset;
            const int col_offset = (pow-1)*nrows;

            double tmp = 0;
#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)
            for(int idx=mat->rowPtr[numa_domain][row]; idx<mat->rowPtr[numa_domain][row+1]; ++idx)
            {
                tmp += mat->val[numa_domain][idx]*x->val[col_offset+mat->col[numa_domain][idx]];
            }
            x->val[(pow)*nrows+row_global] = tmp;
        }
  //  }
}

void plain_spmv_numa(NUMAmat *mat, densemat *x)
{
    ENCODE_TO_VOID_NUMA(mat, NULL, x);
    PLAIN_SPMV_NUMA_KERNEL(0, mat->mat->nrows, 1, voidArg);
    DELETE_ARG();
}

//MKL spmv
void mkl_spmv(densemat* b, sparsemat* mat, densemat* x, bool symm)
{
    int nrows=mat->nrows;
    if(!symm)
    {
        mkl_cspblas_dcsrgemv("N", &nrows, ((double*)mat->val), ((MKL_INT*)mat->rowPtr), ((MKL_INT*) mat->col), ((double*)x->val), ((double*)b->val));
    }
    else
    {
        mkl_cspblas_dcsrsymv("U", &nrows, ((double*)mat->val_symm), ((MKL_INT*)mat->rowPtr_symm), ((MKL_INT*) mat->col_symm), ((double*)x->val), ((double*)b->val));
    }
}

void mkl_spmv(sparsemat* mat, densemat* x, bool symm)
{
    int nrows=mat->nrows;
    double *x_ = x->val;
    double *b_ = &(x->val[mat->nrows]);
    if(!symm)
    {
        mkl_cspblas_dcsrgemv("N", &nrows, ((double*)mat->val), ((MKL_INT*)mat->rowPtr), ((MKL_INT*) mat->col), ((double*)x_), ((double*)b_));
    }
    else
    {
        mkl_cspblas_dcsrsymv("U", &nrows, ((double*)mat->val_symm), ((MKL_INT*)mat->rowPtr_symm), ((MKL_INT*) mat->col_symm), ((double*)x_), ((double*)b_));
    }
}

sparse_matrix_t* mkl_ie_setup(sparsemat* mat, int niter, bool symm)
{
    sparse_matrix_t* A_mkl = new sparse_matrix_t;
    int nrows = mat->nrows;
    struct matrix_descr descr;

    sparse_status_t status;
    if(!symm)
    {
        status = mkl_sparse_d_create_csr(A_mkl, SPARSE_INDEX_BASE_ZERO, nrows, nrows, &(((int*)mat->rowPtr)[0]), &(((int*)mat->rowPtr)[1]), ((int*)mat->col), ((double*)mat->val));
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    }
    else
    {
        mat->computeSymmData();
        status = mkl_sparse_d_create_csr(A_mkl, SPARSE_INDEX_BASE_ZERO, nrows, nrows, &(((int*)mat->rowPtr_symm)[0]), &(((int*)mat->rowPtr_symm)[1]), ((int*)mat->col_symm), ((double*)mat->val_symm));
        descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        descr.mode = SPARSE_FILL_MODE_UPPER;
        descr.diag = SPARSE_DIAG_NON_UNIT;
    }
    if(status!=SPARSE_STATUS_SUCCESS)
    {
        printf("Error ocured in CSR creation of MKL-IE\n");
    }
    status = mkl_sparse_set_mv_hint(*A_mkl, SPARSE_OPERATION_NON_TRANSPOSE, descr, niter);
    if(status!=SPARSE_STATUS_SUCCESS)
    {
        printf("Error ocured in setting hint of MKL-IE\n");
    }
    status = mkl_sparse_optimize(*A_mkl);
    if(status!=SPARSE_STATUS_SUCCESS)
    {
        printf("Error ocured in optimizing phase of MKL-IE\n");
    }
    return A_mkl;
}

void mkl_ie_spmv(densemat* b, sparse_matrix_t* A, densemat* x, bool symm)
{
    struct matrix_descr descr;
    if(!symm)
    {
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    }
    else
    {
        descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        descr.mode = SPARSE_FILL_MODE_UPPER;
        descr.diag = SPARSE_DIAG_NON_UNIT;
    }


    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A, descr, ((double*)x->val), 0, ((double*)b->val));
}

void mkl_ie_spmv(sparse_matrix_t* A, densemat* x, bool symm)
{
    double *x_ = x->val;
    double *b_ = &(x->val[x->nrows]);
    struct matrix_descr descr;
    if(!symm)
    {
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    }
    else
    {
        descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        descr.mode = SPARSE_FILL_MODE_UPPER;
        descr.diag = SPARSE_DIAG_NON_UNIT;
    }

    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A, descr, ((double*)x_), 0, ((double*)b_));
}

void mkl_ie_free(sparse_matrix_t* A)
{
    delete A;
}

inline void MAT_SPMV_KERNEL_only_highest(int start, int end, int pow, int numa_domain, void* args)
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
            const int offset = ((pow-1)%2)*mat->nrows;
#pragma nounroll
#pragma simd vectorlength(VECTOR_LENGTH) reduction(+:tmp)
            for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
            {
                tmp += mat->val[idx]*x->val[offset+mat->col[idx]];
            }
            x->val[(pow%2)*mat->nrows+row] = tmp;
        }
  //  }
}

void matPower_only_highest(sparsemat *A, int power, densemat *x)
{
    RACE::Interface *ce = A->ce;
    ENCODE_TO_VOID(A, NULL, x);
    int race_power_id = ce->registerFunction(&MAT_SPMV_KERNEL_only_highest, voidArg, power);
    {
        ce->executeFunction(race_power_id);
    }
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


