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

    INIT_TIMER(kernel);
    START_TIMER(kernel);

    mat->colorAndPermute(TWO, "RACE", param.cores, param.smt, param.pin);

    STOP_TIMER(kernel);
    double time=GET_TIMER(kernel);
    printf("preprocessing time = %f msec\n", time);
    printf("coloring efficiency = %f\n", mat->colorEff());
    printf("max. stage depth = %d\n", mat->maxStageDepth());

    printf("Finished coloring\n\n");
}
