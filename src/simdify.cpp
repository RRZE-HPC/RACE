#include "simdify.h"

//generate required templates


bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int*rl, int* clp, double* val)
{
    //return testSimdify();
    return simdifyTemplate<double> (simdWidth, C, nrows, col, chunkStart, rl, clp, val);
}

bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int*rl, int* clp, float* val)
{
    return simdifyTemplate<float> (simdWidth, C, nrows, col, chunkStart, rl, clp, val);
}

bool testSimdify()
{
    int C =  4, simdWidth = 4;
    int nrows = 4;

    int *col = new int[16];
    double *val = new double[16];


    for(int i=0; i<nrows; ++i)
    {
        for(int j=0; j<nrows; ++j)
        {
            col[i*nrows+j] = i;
            val[i*nrows+j] = 1;
        }
    }

    int *chunkStart = new int[2];

    chunkStart[0] = 0;
    chunkStart[1] = 17;

    int *rl = new int[nrows];

    for(int i=0; i<nrows; ++i)
    {
        rl[i] = 4;
    }

    int *clp = new int[1];
    clp[0] = 4;

    return simdifyTemplate<double>(simdWidth, C, nrows, col, chunkStart, rl, clp, val);
}
