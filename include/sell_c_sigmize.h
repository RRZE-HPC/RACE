#ifndef RACE_SELL_C_SIGMIZE
#define RACE_SELL_C_SIGMIZE

#include "interface.h"

struct sell_c_sigmize_arg
{
    int simdWidth;
    int C;
    int* col;
    int* chunkStart;
    int* rowLen;
    int* chunkLenPadded;
};


void sell_c_sigmize_Kernel(int start, int end, void *args);
void sell_c_sigmize(int simdWidth, int C, int* col, int* chunkStart, int* rl, int *clp, RACEInterface* ce);

#endif
