#include "kernels.h"
#include <immintrin.h>

#define VECLEN 8
#define NCOLS 64
#include "avx_intrinsics.h"

inline void SPMV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    if(x->ncols == 1)
    {
        if(!mat->isComplex())
        {
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
        else
        {
            for(int row=start; row<end; ++row)
            {
                double tmp_r = 0;
                double tmp_i = 0;
#pragma simd vectorlength(VECLEN) reduction(+:tmp_r,tmp_i)
                for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                {
                    double mat_r = mat->val[2*idx];
                    double mat_i = mat->val[2*idx+1];
                    double x_r = x->val[2*mat->col[idx]];
                    double x_i = x->val[2*mat->col[idx]+1];

                    tmp_r += mat_r*x_r - mat_i*x_i;
                    tmp_i += mat_r*x_i + mat_i*x_r;
                }
                b->val[2*row] = tmp_r;
                b->val[2*row+1] = tmp_i;
            }
        }
    }
    else
    {
        int ncols = x->ncols;
        if(!mat->isComplex())
        {
            double *tmp = new double[ncols];
            for(int row=start; row<end; ++row)
            {
#pragma simd vectorlength(VECLEN)
                for(int col=0; col<ncols; ++col)
                {
                    tmp[col] = 0;
                }
                for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                {
                    double mat_val = mat->val[idx];
                    int col_idx = mat->col[idx];
#pragma simd vectorlength(VECLEN)
                    for(int col=0; col<ncols; ++col)
                    {
                        tmp[col] += mat_val*x->val[col_idx*ncols+col];
                    }
                }
#pragma simd vectorlength(VECLEN)
                 for(int col=0; col<ncols; ++col)
                {
                    b->val[row*ncols+col] = tmp[col];
                }
            }
            delete[] tmp;
        }
        else
        {
            VECREG *tmp = new VECREG[2*ncols/VECLEN];
            for(int row=start; row<end; ++row)
            {
                for(int col=0; col<2*ncols; col+=VECLEN)
                {
                    tmp[col/VECLEN] = SET_ZERO;
                }
#if NCOLS < 8
#pragma unroll(4)
#endif
                for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                {
                    __m128d mat_val_128 = CAST_TO_128(LOAD(&mat->val[2*idx]));
                    VECREG mat_val = SET_FROM_128(mat_val_128, mat_val_128);
                    int col_idx = mat->col[idx];

                    for(int col=0; col<2*ncols; col+=VECLEN)
                    {
                        VECREG x_val = LOAD(&(x->val[(2*(col_idx*ncols))+col]));
                        tmp[col/VECLEN] = ADD(tmp[col/VECLEN], COMPLEX_MUL(mat_val, x_val));
                    }
                }
                for(int col=0; col<2*ncols; col+=VECLEN)
                {
                    STORE(&(b->val[2*(row*ncols)+col]), tmp[col/VECLEN]);
                }
            }
            delete[] tmp;
        }
    }
}

#ifdef __cplusplus
extern "C"
{
#endif
//b=A*x
void spmv_RACE(densemat* b, sparsemat* mat, densemat* x, int iterations)
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
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif
void spmv(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    int start = 0;
    int end = mat->nrows;
    if(x->ncols == 1)
    {
        for(int iter=0; iter<iterations; ++iter)
        {
            if(!mat->isComplex())
            {
#pragma omp parallel for schedule(static)
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
            else
            {
#pragma omp parallel for schedule(static)
                for(int row=start; row<end; ++row)
                {
                    double tmp_r = 0;
                    double tmp_i = 0;
#pragma simd vectorlength(VECLEN) reduction(+:tmp_r,tmp_i)
                    for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                    {
                        double mat_r = mat->val[2*idx];
                        double mat_i = mat->val[2*idx+1];
                        double x_r = x->val[2*mat->col[idx]];
                        double x_i = x->val[2*mat->col[idx]+1];

                        tmp_r += mat_r*x_r - mat_i*x_i;
                        tmp_i += mat_r*x_i + mat_i*x_r;
                    }
                    b->val[2*row] = tmp_r;
                    b->val[2*row+1] = tmp_i;
                }
            }
        }
    }
    else
    {
        const int ncols = NCOLS;//x->ncols;
        for(int iter=0; iter<iterations; ++iter)
        {
            if(!mat->isComplex())
            {
#pragma omp parallel
                {
                    VECREG* tmp = new VECREG[ncols/VECLEN];
#pragma omp for schedule(static)
                    for(int row=start; row<end; ++row)
                    {
                        for(int col=0; col<ncols; col+=VECLEN)
                        {
                            tmp[col/VECLEN] = SET_ZERO;
                        }
#if NCOLS < 16
#pragma unroll(4)
#endif
                        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                        {
                            VECREG mat_val = BROADCAST(&(mat->val[idx]));
                            int col_idx = mat->col[idx];

                            for(int col=0; col<ncols; col+=VECLEN)
                            {
                                VECREG x_val = LOAD(&(x->val[col_idx*ncols+col]));
                                tmp[col/VECLEN] = FMA(mat_val, x_val, tmp[col/VECLEN]);
                            }
                        }
                        for(int col=0; col<ncols; col+=VECLEN)
                        {
                            STORE(&(b->val[row*ncols+col]), tmp[col/VECLEN]);
                        }
                    }
                    delete[] tmp;
                }
            }
            else
            {
                printf("NCOLS = %d\n", NCOLS);
#pragma omp parallel
                {
                    VECREG *tmp = new VECREG[2*ncols/VECLEN];
#pragma omp for schedule(static)
                    for(int row=start; row<end; ++row)
                    {
                        for(int col=0; col<2*ncols; col+=VECLEN)
                        {
                            tmp[col/VECLEN] = SET_ZERO;
                        }
#if NCOLS < 8
#pragma unroll(4)
#endif
                        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                        {
                            __m128d mat_val_128 = CAST_TO_128(LOAD(&mat->val[2*idx]));
                            VECREG mat_val = SET_FROM_128(mat_val_128, mat_val_128);
                            int col_idx = mat->col[idx];

                            for(int col=0; col<2*ncols; col+=VECLEN)
                            {
                                VECREG x_val = LOAD(&(x->val[(2*(col_idx*ncols))+col]));
                                tmp[col/VECLEN] = ADD(tmp[col/VECLEN], COMPLEX_MUL(mat_val, x_val));
                            }
                        }
                        for(int col=0; col<2*ncols; col+=VECLEN)
                        {
                            STORE(&(b->val[2*(row*ncols)+col]), tmp[col/VECLEN]);
                        }
                    }
                    delete[] tmp;
                }
            }
        }
    }
}
#ifdef __cplusplus
}
#endif

inline void SPMTV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    if(x->ncols == 1)
    {
        if(!mat->isComplex())
        {
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
        else
        {
            for(int row=start; row<end; ++row)
            {
                double x_r = x->val[2*row];
                double x_i = x->val[2*row+1];

#pragma simd vectorlength(VECLEN)
                for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                {
                    double mat_r = mat->val[2*idx];
                    double mat_i = -mat->val[2*idx+1]; //conjugate transpose

                    b->val[2*mat->col[idx]] += mat_r*x_r - mat_i*x_i;
                    b->val[2*mat->col[idx]+1] += mat_r*x_i + mat_i*x_r;
                }
            }

        }
    }
    else
    {
        if(!mat->isComplex())
        {
            VECREG* x_row = new VECREG[NCOLS/VECLEN];
            for(int row=start; row<end; ++row)
            {
                for(int col=0; col<NCOLS; col+=VECLEN)
                {
                    x_row[col/VECLEN] = LOAD(&(x->val[row*NCOLS+col]));
                }
#if NCOLS < 8
#pragma unroll(4)
#pragma simd
#endif
                for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                {
                    VECREG mat_val = BROADCAST(&(mat->val[idx]));
                    int col_idx = mat->col[idx];
                    for(int col=0; col<NCOLS; col+=VECLEN)
                    {
                        VECREG b_val = LOAD(&(b->val[col_idx*NCOLS+col]));
                        b_val = FMA(mat_val, x_row[col/VECLEN], b_val);
                        STORE(&(b->val[col_idx*NCOLS+col]), b_val);
                    }
                }
            }
            delete[] x_row;
        }
        else
        {
            VECREG* x_row = new VECREG[2*NCOLS/VECLEN];
            for(int row=start; row<end; ++row)
            {
                for(int col=0; col<2*NCOLS; col+=VECLEN)
                {
                    x_row[col/VECLEN] = LOAD(&(x->val[2*row*NCOLS+col]));
                }
#if NCOLS < 8
#pragma unroll(4)
#endif
                for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
                {
                    __m128d mat_val_128 = CAST_TO_128(LOAD(&mat->val[2*idx]));
                    VECREG mat_val = SET_FROM_128(mat_val_128, mat_val_128);
                    int col_idx = mat->col[idx];

                    for(int col=0; col<2*NCOLS; col+=VECLEN)
                    {
                        VECREG b_val = LOAD(&(b->val[2*col_idx*NCOLS+col]));
                        b_val = ADD(b_val, COMPLEX_MUL_CONJ(mat_val, x_row[col/VECLEN]));
                        STORE(&(b->val[2*col_idx*NCOLS+col]), b_val);
                    }
                }
            }
            delete[] x_row;
        }
    }
}

#ifdef __cplusplus
extern "C"
{
#endif
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
#ifdef __cplusplus
}
#endif

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

#ifdef __cplusplus
extern "C"
{
#endif
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
#ifdef __cplusplus
}
#endif

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

#ifdef __cplusplus
extern "C"
{
#endif
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
#ifdef __cplusplus
}
#endif

inline void SYMM_SPMV_KERNEL(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    if(x->ncols == 1)
    {
        if(!mat->isComplex())
        {
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
        else
        {
            for(int row=start; row<end; ++row)
            {
                double x_row_r = x->val[2*row];
                double x_row_i = x->val[2*row+1];

                double mat_diag_r = mat->val_symm[2*mat->rowPtr_symm[row]];
                double mat_diag_i = mat->val_symm[2*mat->rowPtr_symm[row]+1];

                b->val[2*row] += mat_diag_r*x_row_r - mat_diag_i*x_row_i;
                b->val[2*row+1] += mat_diag_r*x_row_i + mat_diag_i*x_row_r;
                double temp_r = 0;
                double temp_i = 0;

#pragma simd vectorlength(VECLEN) reduction(+:temp_r,temp_i)
                for(int idx=mat->rowPtr_symm[row]+1; idx<mat->rowPtr_symm[row+1]; ++idx)
                {
                    int colIdx = mat->col_symm[idx];
                    double mval_r = mat->val_symm[2*idx];
                    double mval_i = mat->val_symm[2*idx+1];
                    double xval_r = x->val[2*colIdx];
                    double xval_i = x->val[2*colIdx+1];

                    temp_r += mval_r*xval_r - mval_i*xval_i;
                    temp_i += mval_r*xval_i + mval_i*xval_r;
                    //conjgate transpose
                    b->val[2*colIdx] += mval_r*x_row_r + mval_i*x_row_i;
                    b->val[2*colIdx+1] += mval_r*x_row_i - mval_i*x_row_r;
                }
                b->val[2*row]+=temp_r;
                b->val[2*row+1]+=temp_i;
            }
        }
    }
    else
    {
        if(!mat->isComplex())
        {
            VECREG* x_row = new VECREG[NCOLS/VECLEN];
            VECREG* tmp = new VECREG[NCOLS/VECLEN];
            for(int row=start; row<end; ++row)
            {
                for(int col=0; col<NCOLS; col+=VECLEN)
                {
                    x_row[col/VECLEN] = LOAD(&(x->val[row*NCOLS+col]));
                    tmp[col/VECLEN] = SET_ZERO;
                    VECREG b_row = LOAD(&(b->val[row*NCOLS+col]));
                    VECREG m_val = BROADCAST(&(mat->val_symm[mat->rowPtr_symm[row]]));
                    b_row = FMA(m_val, x_row[col/VECLEN], b_row);
                    STORE(&(b->val[row*NCOLS+col]), b_row);
                }
#if NCOLS < 8
#pragma unroll(4)
#endif
                for(int idx=mat->rowPtr_symm[row]+1; idx<mat->rowPtr_symm[row+1]; ++idx)
                {
                    VECREG mval = BROADCAST(&(mat->val_symm[idx]));
                    int colIdx = mat->col_symm[idx];

                    for(int col=0; col<NCOLS; col+=VECLEN)
                    {
                        VECREG x_val = LOAD(&(x->val[colIdx*NCOLS+col]));
                        tmp[col/VECLEN] = FMA(mval, x_val, tmp[col/VECLEN]);
                        VECREG b_val = LOAD(&(b->val[colIdx*NCOLS+col]));
                        b_val = FMA(mval, x_row[col/VECLEN], b_val);
                        STORE(&(b->val[colIdx*NCOLS+col]), b_val);
                    }
                }

                for(int col=0; col<NCOLS; col+=VECLEN)
                {
                    VECREG b_row = LOAD(&(b->val[row*NCOLS+col]));
                    b_row = ADD(b_row, tmp[col/VECLEN]);
                    STORE(&(b->val[row*NCOLS+col]), b_row);
                }
            }
            delete[] x_row;
            delete[] tmp;
        }
        else
        {
            VECREG* x_row = new VECREG[2*NCOLS/VECLEN];
            VECREG* tmp = new VECREG[2*NCOLS/VECLEN];
            for(int row=start; row<end; ++row)
            {
                for(int col=0; col<2*NCOLS; col+=VECLEN)
                {
                    x_row[col/VECLEN] = LOAD(&(x->val[2*row*NCOLS+col]));
                    tmp[col/VECLEN] = SET_ZERO;
                    __m128d mat_val_128 = CAST_TO_128(LOAD(&mat->val_symm[2*mat->rowPtr_symm[row]]));
                    VECREG m_val = SET_FROM_128(mat_val_128, mat_val_128);
                    VECREG bval = LOAD(&(b->val[2*row*NCOLS+col]));
                    bval = ADD(bval, COMPLEX_MUL(m_val,x_row[col/VECLEN]));
                    STORE(&(b->val[2*row*NCOLS+col]), bval);
                }

#if NCOLS < 8
#pragma unroll(4)
#endif
                for(int idx=mat->rowPtr_symm[row]+1; idx<mat->rowPtr_symm[row+1]; ++idx)
                {
                    __m128d mat_val_128 = CAST_TO_128(LOAD(&mat->val_symm[2*idx]));
                    VECREG m_val = SET_FROM_128(mat_val_128, mat_val_128);
                    int colIdx = mat->col_symm[idx];
                    for(int col=0; col<2*NCOLS; col+=VECLEN)
                    {
                        VECREG x_val = LOAD(&(x->val[2*colIdx*NCOLS+col]));
                        tmp[col/VECLEN] = ADD(tmp[col/VECLEN], COMPLEX_MUL(m_val, x_val));

                        VECREG b_val = LOAD(&(b->val[2*colIdx*NCOLS+col]));
                        b_val = ADD(b_val, COMPLEX_MUL_CONJ(m_val, x_row[col/VECLEN]));

                        STORE(&(b->val[2*colIdx*NCOLS+col]), b_val);
                    }
                }

                for(int col=0; col<2*NCOLS; col+=VECLEN)
                {
                    VECREG b_row = LOAD(&(b->val[2*row*NCOLS+col]));
                    b_row = ADD(b_row, tmp[col/VECLEN]);
                    STORE(&(b->val[2*row*NCOLS+col]), b_row);
                }
            }
        }
    }
}

#ifdef __cplusplus
extern "C"
{
#endif
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
#ifdef __cplusplus
}
#endif
