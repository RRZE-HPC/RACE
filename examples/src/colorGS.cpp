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
#include <math.h>

void capitalize(char* beg)
{
    int i = 0;
    while(beg[i])
    {
        beg[i] = toupper(beg[i]);
        ++i;
    }
}

double scaleFn(int i, double scale)
{
    return i*scale;
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
                double errNorm = 0;\
                ERR_CHECK;\
                if(trial==0)\
                {\
                    convergence_history.push_back(std::pair<int, double>(actualIter, errNorm));\
                }\
                if(param.validate) /*Validate is used as verbose, here*/\
                {\
                    double resNorm = 0;\
                    RES_CHECK;\
                    resConvergence_history.push_back(resNorm);\
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

    dist distance=ONE;
    printf("Coloring matrix\n");
    mat->colorAndPermute(distance, std::string(param.colorType), param.cores, param.smt, param.pin);
    printf("Finished coloring\n\n");

    //required for GS kernel
    //mat->makeDiagFirst(0.0, true);
    mat->makeDiagFirst();


    STOP_TIMER(pre_process);
    double pre_process_time = GET_TIMER(pre_process);
    printf("Total pre-processing time = %f s\n", pre_process_time);


    int NROWS = mat->nrows;
    int randInit = false;
    double initVal = 1; ///(double)NROWS;

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
        //std::function<double(int)> initFn = std::bind(&scaleFn, std::placeholders::_1, initVal);
        //xSoln->setFn(initFn);
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
    std::vector<double> resConvergence_history;

    PERF_RUN(color_gs, 2, gs(b, mat, x););

    int lenIter = convergence_history.size();
    bool isDiverging = false;
    //check if convergence or divergence is happening
    if(!isfinite(convergence_history[lenIter-1].second))
    {
        isDiverging=true;
    }
    if((lenIter > 1) && (convergence_history[lenIter-1].second > convergence_history[lenIter-2].second))
    {
        isDiverging=true;
    }

    //if diverging repeat with modified diagonals
    if(isDiverging)
    {
        printf("Given system diverges. Running on a system with modified diagonal values.\n");
        mat->makeDiagFirst(0.0, true);

        b->setVal(0);
        //set b once for LSE, with actual xSoln
        plain_spmv(b, mat, xSoln);

        x->setVal(0);
        actualIter = 0;
        //This macro times and reports performance by running the solver multiple
        //times
        convergence_history.clear();
        PERF_RUN(color_gs, 2, gs(b, mat, x););
    }

    if(param.validate)
    {
        int lenIter = convergence_history.size();
        if(param.convFile == NULL)
        {
            printf("Error history\n");
            printf("#Iter, abs. error, res. error\n");

            for(int i=0; i<lenIter; ++i)
            {
                printf("%d, %.16f, %.16f\n", convergence_history[i].first, convergence_history[i].second, resConvergence_history[i]);
            }
        }
        else
        {
            FILE * fp;
            fp = fopen (param.convFile, "w+");

            fprintf(fp, "#Iter, abs. error, res. error\n");
            for(int i=0; i<lenIter; ++i)
            {
                fprintf(fp, "%d, %.16f, %.16f\n", convergence_history[i].first, convergence_history[i].second, resConvergence_history[i]);
            }
            fclose(fp);
        }
    }

    double resNorm;
    double errNorm;

    RES_CHECK;
    ERR_CHECK;

    printf("Convergence results: resNorm = %e, errNorm = %e, converged iter = %d\n", resNorm, errNorm, actualIter);
    delete mat;
    delete x;
    delete b;
    delete xSoln;
}
