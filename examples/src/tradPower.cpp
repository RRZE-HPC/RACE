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

int main(const int argc, char * argv[])
{
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
    int power = param.iter;
    printf("power = %d\n", power);

    INIT_TIMER(pre_process);
    START_TIMER(pre_process);
    if(param.RCM_flag)
    {
        mat->doRCMPermute(); //Permute();
    }
    STOP_TIMER(pre_process);
    double pre_process_time = GET_TIMER(pre_process);
    printf("Pre-processing time = %f s\n", pre_process_time);


    //mat->writeFile("after_RCM.mtx");
    //mat->prepareForPower(power, param.nodes, param.cache_size*1024*1024, param.cores, param.smt, param.pin);
    //mat->numaInit();

    INFO_PRINT("Matrix statistics");
    INFO_PRINT("Nrows = %d, NNZ = %d, NNZR = %f\n", NROWS, mat->nnz, mat->nnz/((double)NROWS));

    densemat *xTRAD;
    xTRAD=new densemat(NROWS,power+1);
    densemat* xTRAD_0 = xTRAD->view(0,0);

    double initVal = 1/(double)NROWS;
    xTRAD_0->setVal(initVal);

    //determine iterations
    INIT_TIMER(matPower_init);
    START_TIMER(matPower_init);
    for(int iter=0; iter<10; ++iter)
    {
        for(int pow=0; pow<power; ++pow)
        {
            plain_spmv(mat, xTRAD);
        }
    }
    STOP_TIMER(matPower_init);
    double initTime = GET_TIMER(matPower_init);
    int iterations = (int) (1.2*10/initTime);
    //int iterations = 1; //for correctness checking
    printf("Num iterations =  %d\n", iterations);


    xTRAD_0->setVal(initVal);
    //xTRAD->setRand();

    INIT_TIMER(spmv);

#ifdef LIKWID_PERFMON
#pragma omp parallel
        {
            LIKWID_MARKER_START("Trad_power");
        }
#endif
 
    START_TIMER(spmv);
    for(int iter=0; iter<iterations; ++iter)
    {
        for(int pow=0; pow<power; ++pow)
        {
            densemat *x = xTRAD->view(pow,pow);
            plain_spmv(mat,x);
        }
    }
    STOP_TIMER(spmv);
#ifdef LIKWID_PERFMON
#pragma omp parallel
        {
            LIKWID_MARKER_STOP("Trad_power");
        }
#endif
 
    double spmvPowerTime = GET_TIMER(spmv);
    double flops = 2.0*iterations*power*(double)mat->nnz*1e-9;

    printf("SpMV perf. = %f GFlop/s, time = %f s\n", flops/spmvPowerTime, spmvPowerTime);

    delete mat;
    delete xTRAD;

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif

}
