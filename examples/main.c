#include <stdio.h>
#include <omp.h>
#include "mmio.h"
#include "time.h"

#ifdef LIKWID_MEASURE
#include <likwid.h>
#endif
#include "parse.h"
#include "sparsemat.h"
#include "densemat.h"
#include "kernels.h"
#include "timer.h"

int main(const int argc, char * argv[])
{

    int err;
    parser param;
    if(!param.parse_arg(argc, argv))
    {
        printf("Error in reading parameters\n");
    }
    sparsemat* mat = new sparsemat;

    printf("Reading matrix file\n");
    if(!mat->readFile(param.mat_file))
    {
        printf("Error in reading sparse matrix file\n");
    }
    printf("Finished matrix reading\n");

    mat->colorAndPermute(TWO, param.cores, param.smt, param.pin);

    printf("Finished coloring\n");

    int NROWS = mat->nrows;

    densemat *x, *b;
    x=new densemat(NROWS);
    b=new densemat(NROWS);

    x->setRand();
    b->setRand();

    int iter = 100;

    double time = 0;
    double nnz_update = ((double)mat->nnz)*iter*1e-9;

    sleep(1);

    INIT_TIMER(spmv);

    START_TIMER(spmv);
    spmv(b, mat, x, iter);
    STOP_TIMER(spmv);

    time = GET_TIMER(spmv);
    printf("SPMV : %f GFlop/s ; Time = %f\n", 2.0*nnz_update/(time), time);

    sleep(1);

    sleep(1);

    INIT_TIMER(spmtv);

    START_TIMER(spmtv);
    spmtv(b, mat, x, iter);
    STOP_TIMER(spmtv);

    time = GET_TIMER(spmtv);
    printf("SPMTV : %f GFlop/s ; Time = %f\n", 2.0*nnz_update/(time), time);

    sleep(1);

}
