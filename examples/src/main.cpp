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
    printf("Coloring matrix\n");

    mat->colorAndPermute(TWO, param.cores, param.smt, param.pin);

    printf("Finished coloring\n\n");

    int NROWS = mat->nrows;

    densemat *x, *b;
    x=new densemat(NROWS);
    b=new densemat(NROWS);

    x->setRand();
    b->setRand();

    //This macro times and reports performance
    PERF_RUN(spmv,2);
    PERF_RUN(spmtv,2);
    PERF_RUN(kacz,4);
    //diag entry first required for GS
    mat->makeDiagFirst();
    PERF_RUN(gs,2);

    mat->computeSymmData();
    PERF_RUN(symm_spmv,2);


    if(param.validate)
    {
        printf("\n");
        densemat* bSPMV;
        bSPMV = new densemat(NROWS);
        bSPMV->setVal(0);
        b->setVal(0);
        x->setRand();
        //Assuming symmetric matrix provided.
        //Do one SPMV and SPMTV; compare results
        spmv(bSPMV, mat, x, 1);
        spmtv(b, mat, x, 1);

        bool spmtv_flag = checkEqual(bSPMV,b, param.tol);
        if(!spmtv_flag)
        {
            printf("SPMTV failed\n");
        }

        //check symmetric variant also
        //Do one SPMV and symmSPMV; compare results
        b->setVal(0);
        mat->computeSymmData();
        symm_spmv(b, mat, x, 1);
        bool symm_spmv_flag = checkEqual(bSPMV,b, param.tol);
        if(!symm_spmv_flag)
        {
            printf("SYMM SPMV failed\n");
        }

        if(spmtv_flag &&  symm_spmv_flag)
        {
            printf("Validated coloring\n");
        }
        else
        {
            printf("\nValidation failed: Possible reasons:\n");
            printf("1. Tolerance is too low, use option -T to adjust\n");
            printf("2. Matrix not symmetric (right now validation works only for symmetric matrices)\n");
            printf("3. Coloring failed\n\n");
        }
   }
}
