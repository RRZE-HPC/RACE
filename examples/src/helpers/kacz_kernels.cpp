#include "kacz_kernels.h"

inline void KACZ_KERNEL_fission(int start, int end, void* args)
{
    DECODE_FROM_VOID(args);

    for(int row=start; row<end; ++row)
    {
        double rowNorm = 0.0;
        double scale = 0.0;
        int rowLen = mat->rowPtr[row+1]-mat->rowPtr[row];
        int idx = mat->rowPtr[row];
        scale = b->val[row];

        for(int j=0; j<rowLen; ++j)
        {
            scale -= mat->val[idx+j] * x->val[mat->col[idx+j]];
            rowNorm += mat->val[idx+j]*mat->val[idx+j];
        }
        scale /= rowNorm; //omega considered 1

#pragma simd
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            x->val[mat->col[idx]] += scale*mat->val[idx];
        }
    }
}

//Solve for x : A*x=b
void kacz_fission(densemat* b, sparsemat* mat, densemat* x, int iterations)
{
    RACE::Interface *ce = mat->ce;

    ENCODE_TO_VOID(mat,b,x);

    int kaczId = ce->registerFunction(&KACZ_KERNEL_fission, voidArg);

    for(int i=0; i<iterations; ++i)
    {
        ce->executeFunction(kaczId);
    }

    DELETE_ARG();
}

