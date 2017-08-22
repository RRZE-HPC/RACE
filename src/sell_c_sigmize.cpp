#include "sell_c_sigmize.h"

inline void sell_c_sigmize_Kernel(int start, int end, void *args)
{
    sell_c_sigmize_arg *mat = (sell_c_sigmize_arg *)(args);
    int simdWidth = mat->simdWidth;

    int C = mat->C;
    for(int row=start; row<end; ++row)
    {
        int chunk = row/C;
        int rowinchunk = row%C;
        int base_idx = mat->chunkStart[chunk] + rowinchunk;
        for(int j=mat->rowLen[row]; j<mat->chunkLenPadded[chunk]; ++j)
        {
            int idx = base_idx + j*C;
            //row%simdWidth is required; to avoid
            //false sharing within a simd lane;
            //now each 4 row gets different values
            //avoiding such false sharing
            mat->col[idx] = start + (row%simdWidth);
        }
    }
}

void sell_c_sigmize(int simdWidth, int C, int* col, int* chunkStart, int* rl, int *clp, NAMEInterface *ce)
{
    sell_c_sigmize_arg *mat = new sell_c_sigmize_arg;
    mat->simdWidth = simdWidth;
    mat->C = C;
    mat->col = col;
    mat->chunkStart = chunkStart;
    mat->rowLen = rl;
    mat->chunkLenPadded = clp;

    void* voidMat = (void*) (mat);
    int sell_c_sigma_id = ce->registerFunction(&sell_c_sigmize_Kernel, voidMat);
    ce->executeFunction(sell_c_sigma_id);
}
