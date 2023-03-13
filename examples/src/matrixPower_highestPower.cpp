
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
#include <iostream>

//#define VALIDATE_wo_PERM

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

int findMetricId(int group_id, std::string toFind)
{
#ifdef LIKWID_PERFMON
    int numMetrics = perfmon_getNumberOfMetrics(group_id);

    int dataVol_metric_id = numMetrics-1;
    for(int i=0; i<numMetrics; ++i)
    {
        std::string currMetric(perfmon_getMetricName(group_id,i));
        std::cout<<"currMetric = "<<currMetric<<std::endl;
        if(currMetric.find(toFind) != std::string::npos)
        {
            dataVol_metric_id = i;
        }
    }

    return dataVol_metric_id;
#endif
}

void init_likwid()
{
#ifdef LIKWID_PERFMON
    int nthreads;
#pragma omp parallel
    {
#pragma omp single
        nthreads = omp_get_num_threads();
    }

    topology_init();
    CpuTopology_t topo = get_cpuTopology();
    affinity_init();

    int *cpus = (int*)malloc(nthreads * sizeof(int));

    for(int i=0;i<nthreads;i++)
    {
        cpus[i] = topo->threadPool[i].apicId;
    }

    perfmon_init(nthreads, cpus);

    int gid = perfmon_addEventSet("MEM");
    perfmon_setupCounters(gid);
    int mem_metric_id = findMetricId(gid, "Memory data volume\n");
    free(cpus);
#endif
}

void finalize_likwid()
{
#ifdef LIKWID_PERFMON
    //LIKWID_MARKER_CLOSE;
    perfmon_finalize();
    affinity_finalize();
    topology_finalize();
#endif
}

int main(const int argc, char * argv[])
{
//    init_likwid();
#ifdef LIKWID_PERFMON
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

    int NROWS = mat->nrows;
    bool randInit = false;
    double initVal = 1/(double)NROWS;
    densemat *xRAND;
    if(randInit)
    {
        xRAND = new densemat(NROWS);
        xRAND->setRand();
    }
    int power = param.iter;
    printf("power = %d\n", power);

    densemat* xTRAD = NULL;
#ifdef VALIDATE_wo_PERM
    if(param.validate)
    {
        xTRAD=new densemat(NROWS, 2);
        densemat* xTRAD_0 = xTRAD->view(0,0);

        if(randInit)
        {
            xTRAD_0->copyVal(xRAND);
        }
        else
        {
            xTRAD_0->setVal(initVal);
        }

        //now calculate xTRAD in traditional way
        for(int pow=0; pow<power; ++pow)
        {
            //densemat *x = xTRAD->view(pow,pow);
            plain_spmv_only_highest(mat, xTRAD, pow+1);
        }

    }
#endif

    printf("Preparing matrix for power calculation\n");

    //mat->writeFile("matrixBeforeProcessing.mtx");
/*
#ifdef LIKWID_PERFMON
#pragma omp parallel
    {
        LIKWID_MARKER_START("pre_process");
    }
#endif*/
    INIT_TIMER(pre_process);
    START_TIMER(pre_process);
    if(param.RCM_flag)
    {
        mat->doRCM();
    }
    mat->prepareForPower(power, param.cache_size, param.cores, param.smt, param.pin);
    STOP_TIMER(pre_process);
    /*printf("perm = \n");
    for(int i=0; i<NROWS; ++i)
    {
        printf("%d ", mat->finalPerm[i]);
    }
    printf("\n");
    printf("invPerm\n");
    for(int i=0; i<NROWS; ++i)
    {
        printf("%d ", mat->finalInvPerm[i]);
    }
    printf("\n");
*/
    double pre_process_time = GET_TIMER(pre_process);
    printf("Total pre-processing time = %f s\n", pre_process_time);

/*#ifdef LIKWID_PERFMON
#pragma omp parallel
    {
        LIKWID_MARKER_STOP("pre_process");
    }
#endif*/

    //mat->writeFile("matrixAfterProcessing.mtx");

    INFO_PRINT("Matrix statistics");
    INFO_PRINT("Nrows = %d, NNZ = %d, NNZR = %f\n", mat->nrows, mat->nnz, mat->nnz/((double)mat->nrows));


    densemat *x, *xExact;

    //x stores value in the form
    //   x[0],   x[1], ....,   x[nrows-1]
    //  Ax[0],  Ax[1], ....,  Ax[nrows-1]
    // AAx[0], AAx[1], ...., AAx[nrows-1]
    densemat* xRACE=new densemat(NROWS, 2);
    densemat* xRACE_0 = xRACE->view(0,0);
    if(randInit)
    {
        xRACE_0->copyVal(xRAND);
    }
    else
    {
        xRACE_0->setVal(initVal);
    }

    printf("calculation started\n");

    //determine iterations
    INIT_TIMER(matPower_init);
    START_TIMER(matPower_init);
    for(int iter=0; iter<10; ++iter)
    {
       matPower_only_highest(mat, power, xRACE);
    }
    STOP_TIMER(matPower_init);
    double initTime = GET_TIMER(matPower_init);
    int iterations = std::max(1, (int) (1.2*10/initTime));
    //int iterations = 1; //for correctness checking
    printf("Num iterations =  %d\n", iterations);

    double flops = 2.0*power*iterations*(double)mat->nnz*1e-9;

    if(param.validate)
    {
        densemat* xTRAD_perf=new densemat(NROWS, 2);
#ifndef VALIDATE_wo_PERM
        xTRAD=xTRAD_perf;
#endif
        densemat* xTRAD_0 = xTRAD_perf->view(0,0);

       if(randInit)
       {
           xTRAD_0->copyVal(xRAND);
       }
       else
       {
           xTRAD_0->setVal(initVal);
       }

       INIT_TIMER(spmvPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
        {
            LIKWID_MARKER_START("SpMV_power");
        }
#endif
        START_TIMER(spmvPower);
        //now calculate xTRAD in traditional way
        for(int iter=0; iter<iterations; ++iter)
        {

            PAUSE_TIMER(spmvPower);
            if(randInit)
            {
                xTRAD_0->copyVal(xRAND);
            }
            else
            {
                xTRAD_0->setVal(initVal);
            }
            START_TIMER(spmvPower);

            for(int pow=0; pow<power; ++pow)
            {
                //densemat *x = xTRAD_perf->view(pow,pow);
                plain_spmv_only_highest(mat, xTRAD_perf, pow+1);
            }
        }
        STOP_TIMER(spmvPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
        {
            LIKWID_MARKER_STOP("SpMV_power");
        }
#endif
        double spmvPowerTime = GET_TIMER(spmvPower);
        printf("SpMV power perf. = %f GFlop/s, time = %f\n", flops/spmvPowerTime, spmvPowerTime);

#ifdef VALIDATE_wo_PERM
        delete xTRAD_perf;
#endif
        //sleep before going to benchmark
        sleep(1);
    }

    xRACE->setVal(0);
    if(randInit)
    {
        xRACE_0->copyVal(xRAND);
    }
    else
    {
        xRACE_0->setVal(initVal);
    }

    INIT_TIMER(matPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
    {
        LIKWID_MARKER_START("RACE_power");
    }
#endif
    START_TIMER(matPower);
    for(int iter=0; iter<iterations; ++iter)
    {
        PAUSE_TIMER(matPower);
        if(randInit)
        {
            xRACE_0->copyVal(xRAND);
        }
        else
        {
            xRACE_0->setVal(initVal);
        }
        START_TIMER(matPower);
        matPower_only_highest(mat, power, xRACE);
    }
    STOP_TIMER(matPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
    {
        LIKWID_MARKER_STOP("RACE_power");
    }
#endif
    double RACEPowerTime = GET_TIMER(matPower);
    printf("RACE power perf. = %f GFlop/s, time = %f s\n", flops/RACEPowerTime, RACEPowerTime);


    if(param.validate)
    {
        //permute the trad vector for comparison
#ifdef VALIDATE_wo_PERM
        densemat* xTRAD_permuted = mat->permute_densemat(xTRAD);
#else
        densemat* xTRAD_permuted = xTRAD;
#endif


        xRACE->setVal(0);
        if(randInit)
        {
            xRACE_0->copyVal(xRAND);
        }
        else
        {
            xRACE_0->setVal(initVal);
        }

#ifdef VALIDATE_wo_PERM
        densemat* xRACE_permuted = mat->permute_densemat(xRACE);
#else
        densemat* xRACE_permuted = xRACE;
#endif

        //only one iterations
        matPower_only_highest(mat, power, xRACE_permuted);
        /*        for(int i=0; i<10; ++i)
        {
            for(int j=0; j<xRACE->ncols; ++j)
            {
                printf("(%.8f, %.8f) ", xRACE->val[j*xRACE->nrows+i], xTRAD_permuted->val[j*xRACE->nrows+i]);
            }
            printf("\n");
        }*/

        bool flag = checkEqual(xTRAD_permuted, xRACE_permuted, param.tol);
        if(!flag)
        {
            printf("Power calculation failed\n");
        }
        else
        {
            printf("Power calculation success\n");
        }
#ifdef VALIDATE_wo_PERM
        delete xTRAD_permuted;
        delete xRACE_permuted;
#endif
        delete xTRAD;
    }

    delete mat;
    delete xRACE;
    if(randInit)
    {
        delete xRAND;
    }
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif
}
