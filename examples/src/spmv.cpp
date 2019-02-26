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

    mat->numaInit();

    int NROWS = mat->nrows;
    int iterations = param.iter;

    densemat *x;
    x=new densemat(NROWS,2);

    x->setRand();

    INIT_TIMER(spmv);
    START_TIMER(spmv);
    for(int iter=0; iter<iterations; ++iter)
    {
        plain_spmv(mat,x);
    }
    STOP_TIMER(spmv);
    double spmvPowerTime = GET_TIMER(spmv);
    double flops = 2.0*iterations*(double)mat->nnz*1e-9;

    printf("SpMV perf. = %f, time = %f\n", flops/spmvPowerTime, spmvPowerTime);

    delete mat;
    delete x;
}
