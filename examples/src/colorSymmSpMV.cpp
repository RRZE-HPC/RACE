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
#include "quantile.hpp"

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

#define PERF_RUN(kernel, flopPerNnz, ...)\
{\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
/*_Pragma("omp parallel")\
    {\
        LIKWID_MARKER_REGISTER(capsKernel);\
    }\*/\
    int num_trials=10;\
    std::vector<double> time(num_trials, 0);\
    double nnz_update = ((double)mat->nnz)*iterations*1e-9;\
    INIT_TIMER(kernel);\
    _Pragma("omp parallel")\
    {\
        LIKWID_MARKER_START(capsKernel);\
    }\
    for(int trial=0; trial<num_trials; ++trial)\
    {\
        START_TIMER(kernel);\
        for(int iter=0; iter<iterations; ++iter)\
        {\
            __VA_ARGS__;\
            /*kernel(b, mat, x, symm);*/\
        }\
        STOP_TIMER(kernel);\
        time[trial] = GET_TIMER(kernel);\
    }\
    _Pragma("omp parallel")\
    {\
        LIKWID_MARKER_STOP(capsKernel);\
    }\
    std::vector<double> time_quantiles = Quantile(time, {0, 0.25, 0.5, 0.75, 1});\
   /*printf("Obtained Perf of %s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);*/\
    printf("Obtained Perf of %s : ", capsKernel);\
    Quantile_print(time_quantiles, flopPerNnz*nnz_update, true);\
    printf(" GFlop/s; Time : ");\
    Quantile_print(time_quantiles);\
    printf(" sec\n");\
    free(capsKernel);\
}\

#else

#define PERF_RUN(kernel, flopPerNnz, ...)\
{\
    int num_trials=10;\
    std::vector<double> time(num_trials, 0);\
    double nnz_update = ((double)mat->nnz)*iterations*1e-9;\
    INIT_TIMER(kernel);\
    for(int trial=0; trial<num_trials; ++trial)\
    {\
        START_TIMER(kernel);\
        for(int iter=0; iter<iterations; ++iter)\
        {\
            __VA_ARGS__;\
            /*kernel(b, mat, x, symm);*/\
        }\
        STOP_TIMER(kernel);\
        time[trial] = GET_TIMER(kernel);\
    }\
    std::vector<double> time_quantiles = Quantile(time, {0, 0.25, 0.5, 0.75, 1});\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    /*printf("Obtained Perf of %s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);*/\
    printf("Obtained Perf of %s : ", capsKernel);\
    Quantile_print(time_quantiles, flopPerNnz*nnz_update, true);\
    printf(" GFlop/s; Time : ");\
    Quantile_print(time_quantiles);\
    printf(" sec\n");\
    free(capsKernel);\
}\

#endif

int main(const int argc, char * argv[])
{
    /*int nthreads=1;
#pragma omp parallel
    {
        nthreads=omp_get_num_threads();
    }
    printf("Running with %d threads\n", nthreads);
*/
    /*int curPID = getpid();
    char *curPID_str;
    asprintf(&curPID_str, "%d", curPID);
    printf("Current PID = %s\n", curPID_str);
    setenv("LIKWID_PERF_PID", curPID_str, 1);
    */
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

    INIT_TIMER(pre_process);
    START_TIMER(pre_process);
    if(param.RCM_flag)
    {
        mat->doRCM();
    }

    dist distance=TWO;

    char *dist_env = (getenv("COLOR_DISTANCE"));
    if(dist_env)
    {
        int dist_numerical = atoi(dist_env);
        if(dist_numerical == 1)
        {
            distance=ONE;
        }
    }
    printf("Coloring matrix\n");
    mat->colorAndPermute(distance, std::string(param.colorType), param.cores, param.smt, param.pin);
    printf("Finished coloring\n\n");

    STOP_TIMER(pre_process);
    double pre_process_time = GET_TIMER(pre_process);
    printf("Total pre-processing time = %f s\n", pre_process_time);


    int NROWS = mat->nrows;
    int NNZ = mat->nnz;
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
       plain_spmv(b, mat, x);
    }
    STOP_TIMER(init_iter);
    double initTime = GET_TIMER(init_iter);
    int iterations = std::max(1, (int) (0.1*2*10/initTime));
    //int iterations = 1; //for correctness checking
    printf("Num iterations =  %d\n", iterations);

    //This macro times and reports performance
    b->setVal(0);
    PERF_RUN(omp_spmv,2, plain_spmv(b, mat, x););
    b->setVal(0);
    PERF_RUN(color_spmv,2, color_spmv(b, mat, x););
    //mat->printTree();

    mat->computeSymmData();

    b->setVal(0);
    printf("NROWS = %d, NNZ = %d, NNZR = %f, NNZ_symm = %d, NNZR_symm = %f\n", NROWS, NNZ, NNZ/((double)NROWS), mat->nnz_symm, mat->nnz_symm/((double)NROWS));
    PERF_RUN(color_symm_spmv,2, symm_spmv(b, mat, x););

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

        //Do one SPMV and symmSPMV; compare results
        plain_spmv(bSPMV, mat, x);

        b->setVal(0);
        color_spmv(b, mat, x);
        printf("Checking max deviation between SpMV and coloring run SpMV\n");
        findMaxDeviations(bSPMV, b);
        bool spmv_flag = checkEqual(bSPMV, b, param.tol, relativeCheck);
        if(!spmv_flag)
        {
            printf("Color SPMV failed\n");
        }
        printf("check done\n");

        b->setVal(0);
        mat->computeSymmData();
        symm_spmv(b, mat, x);
        printf("Checking max deviation between SpMV and coloring run SymmSpMV\n");
        findMaxDeviations(bSPMV, b);
        bool symm_spmv_flag = checkEqual(bSPMV, b, param.tol, relativeCheck);
        if(!symm_spmv_flag)
        {
            printf("Color SYMM SPMV failed\n");
        }
        printf("check done\n");

        if(spmv_flag && symm_spmv_flag)
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

#ifdef LIKWID_MEASURE
    LIKWID_MARKER_CLOSE;
#endif


    delete mat;
    delete x;
    delete b;
    if(xRAND)
    {
        delete xRAND;
    }

}
