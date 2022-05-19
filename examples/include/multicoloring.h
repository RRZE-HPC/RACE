#ifndef _MULTICOLORING_PERM_
#define _MULTICOLORING_PERM_

#include "sparsemat.h"

class multicoloring
{
    int* initPerm;
    int* initInvPerm;

    public:
    sparsemat* mat;
    int colorDist;
    int colorBlockSize_out;
    int ncolors_out;
    int* colorPtr_out;
    int* partPtr_out;
    int* perm_out;
    int* invPerm_out;

    multicoloring(sparsemat* mat, int colorDist, int* initPerm, int* initInvPerm);
    bool doMC();
    bool doABMC();
};

#endif
