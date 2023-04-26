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

#ifdef LIKWID_MEASURE

#define PERF_RUN(kernel, flopPerNnz)\
{\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    int iter = param.iter;\
    double time = 0;\
    double nnz_update = ((double)mat->nnz)*iterations*1e-9;\
    sleep(1);\
    INIT_TIMER(kernel);\
    _Pragma("omp parallel")\
    {\
        LIKWID_MARKER_START(capsKernel);\
    }\
    START_TIMER(kernel);\
    for(int iter=0; iter<iterations; ++iter)\
    {\
        kernel(b, mat, x);\
    }\
    STOP_TIMER(kernel);\
    _Pragma("omp parallel")\
    {\
        LIKWID_MARKER_STOP(capsKernel);\
    }\
    time = GET_TIMER(kernel);\
    printf("%10s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);\
    free(capsKernel);\
}\

#else

#define PERF_RUN(kernel, flopPerNnz)\
{\
    int iter = param.iter;\
    double time = 0;\
    double nnz_update = ((double)mat->nnz)*iterations*1e-9;\
    sleep(1);\
    INIT_TIMER(kernel);\
    START_TIMER(kernel);\
    for(int iter=0; iter<iterations; ++iter)\
    {\
        kernel(b, mat, x);\
    }\
    STOP_TIMER(kernel);\
    time = GET_TIMER(kernel);\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    printf("%10s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);\
    free(capsKernel);\
}\

#endif

int main(const int argc, char * argv[])
{
#ifdef LIKWID_MEASURE
    LIKWID_MARKER_INIT;
#endif

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

    if(param.RCM_flag)
    {
        mat->doRCMPermute();
    }

    //mat->numaInit();
    int NROWS = mat->nrows;
    int randInit = false;
    double initVal = 1/(double)NROWS;
    int NNZ = mat->nnz;
    int NNZ_symm = ((NNZ-NROWS)/2)+NROWS;
    printf("Nrows = %d, NNZ = %d, NNZR = %f, NNZ_symm = %d, NNZR_symm = %f\n", NROWS, NNZ, NNZ/((double)NROWS), NNZ_symm, NNZ_symm/((double)NROWS));
    densemat *x, *b;
    x=new densemat(NROWS);
    b=new densemat(NROWS);

    densemat *xRAND;
    if(randInit)
    {
        xRAND = new densemat(NROWS);
        xRAND->setRand();
    }

    b->setVal(0);
    if(randInit)
    {
        x->copyVal(xRAND);
    }
    else
    {
        x->setVal(initVal);
    }


    //determine iterations
    INIT_TIMER(init_iter);
    START_TIMER(init_iter);
    for(int iter=0; iter<10; ++iter)
    {
       plain_spmv(b, mat, x);
    }
    STOP_TIMER(init_iter);
    double initTime = GET_TIMER(init_iter);
    int iterations = std::max(1, (int) (2*10/initTime));
    //int iterations = 1; //for correctness checking
    printf("Num iterations =  %d\n", iterations);


    //This macro times and reports performance
    PERF_RUN(plain_spmv,2);

#ifdef LIKWID_MEASURE
    LIKWID_MARKER_CLOSE;
#endif

    delete x;
    delete b;
    delete mat;
}
