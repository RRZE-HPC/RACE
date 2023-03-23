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

#define RES_CHECK\
    plain_spmv(Ax, mat, x);\
    res->axpby(b, Ax, 1, -1);\
    resNorm = (res->dot(res))[0];\
    resNorm = sqrt(resNorm);\

#define ERR_CHECK\
    err->axpby(x, xSoln, 1, -1);\
    errNorm = (err->dot(err))[0];\
    errNorm = sqrt(errNorm);\


#define PERF_RUN(kernel, flopPerNnz, ...)\
{\
    std::vector<double> time; /*(num_trials, 0);*/\
    INIT_TIMER(kernel);\
    bool converged = false;\
    actualIter = 0;\
    for(int trial=0; trial<num_trials; ++trial)\
    {\
        actualIter = 0;\
        x->setVal(0); /*restart*/\
        START_TIMER(kernel);\
        INIT_TIMER(convCheck);\
        for(int iter=0; iter<iterations; ++iter)\
        {\
            __VA_ARGS__;\
            ++actualIter;\
            /*find error norm every 10 iter*/\
            START_TIMER(convCheck);\
            if((actualIter%10) == 0)\
            {\
                double errNorm;\
                ERR_CHECK;\
                if(trial==0)\
                {\
                    convergence_history.push_back(std::pair<int, double>(actualIter, errNorm));\
                }\
                if(errNorm < param.convTol)\
                {\
                    converged = true;\
                    break;\
                }\
            }\
            PAUSE_TIMER(convCheck);\
        }\
        /*if(converged)\
        {\
            break;\
        }*/\
        STOP_TIMER(convCheck);\
        STOP_TIMER(kernel);\
        time.push_back(GET_TIMER(kernel)-GET_TIMER(convCheck));\
    }\
    std::vector<double> time_quantiles = Quantile(time, {0, 0.25, 0.5, 0.75, 1});\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    /*printf("Obtained Perf of %s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);*/\
    printf("Obtained Perf of %s : ", capsKernel);\
    double nnz_update = ((double)mat->nnz)*actualIter*1e-9;\
    Quantile_print(time_quantiles, flopPerNnz*nnz_update, true);\
    printf(" GFlop/s; Time : ");\
    Quantile_print(time_quantiles);\
    printf(" sec\n");\
    free(capsKernel);\
}\


int main(const int argc, char * argv[])
{

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
    printf("Coloring matrix\n");
    mat->colorAndPermute(distance, std::string(param.colorType), param.cores, param.smt, param.pin);
    printf("Finished coloring\n\n");

    STOP_TIMER(pre_process);
    double pre_process_time = GET_TIMER(pre_process);
    printf("Total pre-processing time = %f s\n", pre_process_time);


    int NROWS = mat->nrows;
    int randInit = false;
    double initVal = 1/(double)NROWS;

    densemat *x, *b;
    x=new densemat(NROWS);
    b=new densemat(NROWS);

    densemat *xSoln=new densemat(NROWS);
    if(randInit)
    {
        xSoln->setRand();
    }
    else
    {
        xSoln->setVal(initVal);
    }

    int iterations;
    if(param.iter == -1)
    {
        //determine iterations
        INIT_TIMER(init_iter);
        START_TIMER(init_iter);
        for(int iter=0; iter<10; ++iter)
        {
            plain_spmv(b, mat, x);
        }
        STOP_TIMER(init_iter);
        double initTime = GET_TIMER(init_iter);
        iterations = std::max(20, (int) (0.1*2*10/initTime));
    }
    else
    {
        iterations = param.iter;
    }
    int num_trials=10;
    printf("Num iterations =  %d\n", iterations);

    b->setVal(0);
    //set b once for LSE, with actual xSoln
    plain_spmv(b, mat, xSoln);

    //set convergence check vectors;
    densemat *Ax = new densemat(NROWS);
    densemat *res = new densemat(NROWS);
    densemat *err = new densemat(NROWS);

    x->setVal(0);
    int actualIter = 0;
    //This macro times and reports performance by running the solver multiple
    //times
    std::vector<std::pair<int, double>> convergence_history;
    PERF_RUN(color_kacz, 4, kacz(b, mat, x););

    if(param.validate)
    {
        int lenIter = convergence_history.size();
        if(param.convFile == NULL)
        {
            printf("Error history\n");
            printf("#Iter, abs. error\n");

            for(int i=0; i<lenIter; ++i)
            {
                printf("%d, %.16f\n", convergence_history[i].first, convergence_history[i].second);
            }
        }
        else
        {
            FILE * fp;
            fp = fopen (param.convFile, "w+");

            fprintf(fp, "#Iter, abs. error\n");
            for(int i=0; i<lenIter; ++i)
            {
                fprintf(fp, "%d, %.16f\n", convergence_history[i].first, convergence_history[i].second);
            }
            fclose(fp);

        }
    }

    double resNorm;
    RES_CHECK;
    double errNorm;
    ERR_CHECK;

    printf("Convergence results: resNorm = %e, errNorm = %e, converged iter = %d\n", resNorm, errNorm, actualIter);
    delete mat;
    delete x;
    delete b;
    delete xSoln;
}
