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

#define PERF_RUN(kernel, flopPerNnz, ...)\
{\
    double time = 0;\
    double nnz_update = ((double)mat->nnz)*iterations*1e-9;\
    INIT_TIMER(kernel);\
    START_TIMER(kernel);\
    for(int iter=0; iter<iterations; ++iter)\
    {\
        __VA_ARGS__;\
        /*kernel(b, mat, x, symm);*/\
    }\
    STOP_TIMER(kernel);\
    time = GET_TIMER(kernel);\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    printf("%s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);\
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

    INIT_TIMER(pre_process);
    START_TIMER(pre_process);
    if(param.RCM_flag)
    {
        mat->doRCMPermute();
    }

    STOP_TIMER(pre_process);
    double pre_process_time = GET_TIMER(pre_process);
    printf("Total pre-processing time = %f s\n", pre_process_time);


    int NROWS = mat->nrows;
    int randInit = false;
    double initVal = 1/(double)NROWS;

    densemat *x, *b;
    x=new densemat(NROWS);
    b=new densemat(NROWS);

    densemat *xRAND=NULL;
    if(randInit)
    {
        xRAND = new densemat(NROWS);
        xRAND->setRand();
    }

    if(randInit)
    {
        x->copyVal(xRAND);
    }
    else
    {
        x->setVal(initVal);
    }

    b->setVal(0);
    //determine iterations
    INIT_TIMER(init_iter);
    START_TIMER(init_iter);
    for(int iter=0; iter<10; ++iter)
    {
       mkl_spmv(b, mat, x);
    }
    STOP_TIMER(init_iter);
    double initTime = GET_TIMER(init_iter);
    int iterations = std::max(1, (int) (2*10/initTime));
    //int iterations = 1; //for correctness checking
    printf("Num iterations =  %d\n", iterations);

    //This macro times and reports performance
    b->setVal(0);
    PERF_RUN(mkl_spmv,2, mkl_spmv(b, mat, x););
    mat->computeSymmData();
    b->setVal(0);
    PERF_RUN(mkl_symm_spmv,2, mkl_spmv(b, mat, x, true););

    //now run mkl-IE
    sparse_matrix_t* A_mkl=mkl_ie_setup(mat, iterations);
    b->setVal(0);
    PERF_RUN(mkl_ie_spmv, 2, mkl_ie_spmv(b, A_mkl, x););
    sparse_matrix_t* A_symm_mkl=mkl_ie_setup(mat, iterations, true);
    b->setVal(0);
    PERF_RUN(mkl_ie_symm_spmv, 2, mkl_ie_spmv(b, A_symm_mkl, x, true););

    if(param.validate)
    {
        printf("\n");
        densemat* bSPMV;
        bSPMV = new densemat(NROWS);
        bSPMV->setVal(0);
        if(randInit)
        {
            x->copyVal(xRAND);
        }
        else
        {
            x->setVal(initVal);
        }

        bool relativeCheck = false;
        if(randInit)
        {
            relativeCheck = true; //for random normally you got very big numbers, so check relative error
        }
        //check symmetric variant also
        //Do one SPMV and symmSPMV; compare results
        mkl_spmv(bSPMV, mat, x, false);
        mat->computeSymmData();
        b->setVal(0);
        mkl_spmv(b, mat, x, true);
        bool symm_spmv_flag = checkEqual(bSPMV, b, param.tol, relativeCheck);
        if(!symm_spmv_flag)
        {
            printf("SYMM SPMV failed\n");
        }
        b->setVal(0);
        mkl_ie_spmv(b, A_symm_mkl, x, true);
        bool ie_symm_spmv_flag = checkEqual(bSPMV, b, param.tol, relativeCheck);
        if(!ie_symm_spmv_flag)
        {
            printf("SYMM SPMV IE failed\n");
        }

 
        if(symm_spmv_flag && ie_symm_spmv_flag)
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
        delete bSPMV;
   }

    delete mat;
    delete x;
    delete b;
    if(xRAND)
    {
        delete xRAND;
    }
    mkl_ie_free(A_mkl);
    mkl_ie_free(A_symm_mkl);

}
