#ifndef RACE_KACZ_KERNELS_H
#define RACE_KACZ_KERNELS_H

#include "sparsemat.h"
#include "densemat.h"
#include "kernels.h"

void kacz_fission(densemat* b, sparsemat* mat, densemat* x, int iter);

#endif
