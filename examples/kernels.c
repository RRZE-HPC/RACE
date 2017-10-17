#include "kernels.h"

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
    RACEInterface *ce = mat->ce;

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
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            b->val[mat->col[idx]] += mat->val[idx]*x->val[row];
        }
    }
}


//b=A'*x
void spmtv(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    RACEInterface *ce = mat->ce;

    ENCODE_TO_VOID(mat,b,x);

    int spmtvId = ce->registerFunction(&SPMTV_KERNEL, voidArg);

    for(int i=0; i<iterations; ++i)
    {
        ce->executeFunction(spmtvId);
    }

    DELETE_ARG();
}
