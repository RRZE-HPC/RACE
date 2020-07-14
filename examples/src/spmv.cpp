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


void capitalize(char* beg)
{
    int i = 0;
    while(beg[i])
    {
        beg[i] = toupper(beg[i]);
        ++i;
    }
}

#define PERF_RUN(kernel, flopPerNnz)\
{\
    int iter = param.iter;\
    double time = 0;\
    double nnz_update = ((double)mat->nnz)*iter*1e-9;\
    sleep(1);\
    INIT_TIMER(kernel);\
    START_TIMER(kernel);\
    kernel(b, mat, x, iter);\
    STOP_TIMER(kernel);\
    time = GET_TIMER(kernel);\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    printf("%10s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);\
    free(capsKernel);\
}\

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

    mat->doRCMPermute();
    //mat->numaInit();
    int NROWS = mat->nrows;

    densemat *x, *b;
    x=new densemat(NROWS);
    b=new densemat(NROWS);

    x->setRand();
    b->setRand();

    //This macro times and reports performance
    PERF_RUN(plain_spmv,2);

    delete x;
    delete b;
    delete mat;
}
