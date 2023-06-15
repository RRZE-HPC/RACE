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

    //At this point matrix is created and all pre-processing is done

    int NROWS = mat->nrows;
    int randInit = false;
    double initVal = 1/(double)NROWS;

    densemat *x, *b;
    x=new densemat(NROWS);
    b=new densemat(NROWS);

    //store solution vector for testing
    densemat *xSoln=new densemat(NROWS);
    if(randInit)
    {
        xSoln->setRand();
    }
    else
    {
        xSoln->setVal(initVal);
    }

    int iterations = param.iter;
    int num_trials=10;
    printf("Num iterations =  %d\n", iterations);

    b->setVal(0);
    //set b once for LSE, with actual xSoln
    plain_spmv(b, mat, xSoln);

    //set convergence check vectors;
    densemat *Ax = new densemat(NROWS);
    densemat *res = new densemat(NROWS);
    densemat *err = new densemat(NROWS);

    //This is the unknown solution vector x
    //The vector should converge to xSoln
    x->setVal(0);

    //write your CGMN algorithm here
    //Vectors required for the algorithm can be created with the class densemat (See the definition in examples/helpers/densemat.cpp)
    //Also the operation on vector can be seen in the same class.
    //
    //Sparsematrix A is already available in the variable mat. Operations on
    //sparsemat can be seen in examples/helpers/kernels.cpp
    int actualIter = 0; //to store the iterations required to converge



    //Finally, calculate the residual and error norma and report
    double resNorm;
    RES_CHECK;
    double errNorm;
    ERR_CHECK;

    printf("Convergence results: resNorm = %.8f, errNorm = %.8f, converged iter = %d\n", resNorm, errNorm, actualIter);
    delete mat;
    delete x;
    delete b;
    delete xSoln;
}
