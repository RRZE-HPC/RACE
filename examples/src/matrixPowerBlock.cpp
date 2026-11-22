//BCSR matrix power A^p x with RACE level blocking on the node graph.
//Input must be node-ordered: node i = point rows b*i .. b*i+b-1. Flow: read point CRS -> convertToBlockNodes(b)
//-> optional RCM (-R) on the node graph -> RACE with the cache size scaled to the real bytes per block entry.
//Same options and output lines as mtxPower. GFlop/s count 2 flops per stored point entry (explicit zeros of the
//blocks not counted), so they compare directly with mtxPower on the same file. "SpMV power perf." = plain BCSR
//SpMV, "RACE power perf." = BCSR + RACE, so the RACE/SpMV ratio of this binary is BCSR + RACE vs BCSR SpMV.
//Validation compares A^p x against point CRS SpMV in the original numbering.
//Environment: MTXPOWER_BLOCK_SIZE (default 3), MTXPOWER_ITERS (as mtxPower).
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <unistd.h>

#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif
#include "parse.h"
#include "sparsemat.h"
#include "densemat.h"
#include "kernels.h"
#include "timer.h"
#include <algorithm>

int main(const int argc, char * argv[])
{
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_INIT;
#endif
    parser param;
    if(!param.parse_arg(argc, argv))
    {
        printf("Error in reading parameters\n");
    }
    int b = 3;
    const char* bsEnv = getenv("MTXPOWER_BLOCK_SIZE");
    if(bsEnv && atoi(bsEnv) > 0)
    {
        b = atoi(bsEnv);
    }
    printf("block size = %d\n", b);

    sparsemat* mat = new sparsemat;
    printf("Reading matrix file\n");
    if(!mat->readFile(param.mat_file))
    {
        printf("Error in reading sparse matrix file\n");
        return 1;
    }
    const int NROWS = mat->nrows;   //point rows
    const double initVal = 1/(double)NROWS;
    int power = param.iter;
    printf("power = %d\n", power);

    //point copy in the original numbering, for the validation only
    sparsemat* pmat = NULL;
    if(param.validate)
    {
        pmat = new sparsemat;
        pmat->basicDeepCopy(mat);
    }

    printf("Preparing matrix for power calculation\n");
    INIT_TIMER(pre_process);
    START_TIMER(pre_process);
    if(!mat->convertToBlockNodes(b))
    {
        printf("Error in BCSR conversion\n");
        return 1;
    }
    if(param.RCM_flag)
    {
        mat->doRCM();   //on the node graph (block pattern); prepareForPower passes it to RACE as initPerm
    }

    //RACE sizes level groups with bytePerNNZ = (1+p/nnzr)*8 + (1+1/nnzr)*4 (one value and one column index per
    //entry, p vector values per row; mtxPower::getBytePerNNZ). A block entry costs b*b*8 + 4 B and a node row
    //holds b*p vector values, so the cache is scaled by model/true bytes per entry.
    const double nnzr = mat->nnz/(double)mat->nrows;
    const double modelBytes = (1.0 + power/nnzr)*sizeof(double) + (1.0 + 1.0/nnzr)*sizeof(int);
    const double trueBytes = ((double)b*b + (double)power*b/nnzr)*sizeof(double) + (1.0 + 1.0/nnzr)*sizeof(int);
    const double cacheScale = modelBytes/trueBytes;
    const double cacheEff = param.cache_size*cacheScale;
    printf("BCSR cache model: blocks per node row = %f, RACE bytes/entry = %f, true bytes/entry = %f, "
           "cache %f MB -> passed to RACE %f MB (scale %f)\n",
           nnzr, modelBytes, trueBytes, param.cache_size, cacheEff, cacheScale);

    mat->prepareForPower(power, cacheEff, param.cores, param.smt, param.pin);
    STOP_TIMER(pre_process);
    printf("Total pre-processing time = %f s\n", GET_TIMER(pre_process));

    INFO_PRINT("Matrix statistics");
    INFO_PRINT("Nrows = %d, NNZ = %lld, NNZR = %f (point); node rows = %d, blocks = %d\n",
               NROWS, mat->pointNnz, mat->pointNnz/((double)NROWS), mat->nrows, mat->nnz);

    densemat* xRACE = new densemat(NROWS, power+1);
    densemat* xRACE_0 = xRACE->view(0,0);
    xRACE_0->setVal(initVal);

    printf("calculation started\n");

    //determine iterations (as mtxPower)
    INIT_TIMER(matPower_init);
    START_TIMER(matPower_init);
    for(int iter=0; iter<10; ++iter)
    {
        matPowerBlock(mat, power, xRACE);
    }
    STOP_TIMER(matPower_init);
    double initTime = GET_TIMER(matPower_init);
    int iterations = std::max(1, (int) (1.2*10/initTime));
    const char* fixedIters = getenv("MTXPOWER_ITERS");
    if(fixedIters && atoi(fixedIters) > 0)
    {
        iterations = atoi(fixedIters);
    }
    printf("Num iterations =  %d\n", iterations);

    //useful flops: 2 per stored point entry (explicit zeros of the blocks not counted)
    const double flops = 2.0*power*iterations*(double)mat->pointNnz*1e-9;
    printf("BCSR flops incl. explicit zeros = %f x useful\n", (double)b*b*mat->nnz/(double)mat->pointNnz);

    if(param.validate)
    {
        //baseline: plain BCSR SpMV, one power after the other, same (RACE-permuted) block matrix
        densemat* xTRAD_perf = new densemat(NROWS, power+1);
        densemat* xTRAD_0 = xTRAD_perf->view(0,0);
        xTRAD_0->setVal(initVal);
        INIT_TIMER(spmvPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
        {
            LIKWID_MARKER_START("SpMV_power");
        }
#endif
        START_TIMER(spmvPower);
        for(int iter=0; iter<iterations; ++iter)
        {
            for(int pow=0; pow<power; ++pow)
            {
                densemat *x = xTRAD_perf->view(pow,pow);
                plain_spmv_block(mat, x);
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
        delete xTRAD_perf;
        sleep(1);
    }

    xRACE->setVal(0);
    xRACE_0->setVal(initVal);

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
        matPowerBlock(mat, power, xRACE);
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
        //reference: point CRS, original numbering
        densemat* xRef = new densemat(NROWS, power+1);
        xRef->view(0,0)->setVal(initVal);
        for(int pow=0; pow<power; ++pow)
        {
            plain_spmv(pmat, xRef->view(pow,pow));
        }

        //one BCSR + RACE power, mapped back: new node i holds old node finalPerm[i]
        xRACE->setVal(0);
        xRACE_0->setVal(initVal);
        matPowerBlock(mat, power, xRACE);
        densemat* xMapped = new densemat(NROWS, power+1);
        const int nnodes = mat->nrows;
#pragma omp parallel for schedule(static)
        for(int i=0; i<nnodes; ++i)
        {
            int old = mat->finalPerm[i];
            for(int pow=0; pow<=power; ++pow)
            {
                for(int k=0; k<b; ++k)
                {
                    xMapped->val[(size_t)pow*NROWS + (size_t)b*old + k] = xRACE->val[(size_t)pow*NROWS + (size_t)b*i + k];
                }
            }
        }

        double maxAbs = 0, maxRel = 0;
        long long nExact = 0, nTotal = (long long)(power+1)*NROWS;
        for(long long j=0; j<nTotal; ++j)
        {
            double d = fabs(xRef->val[j] - xMapped->val[j]);
            if(d == 0.0)
            {
                ++nExact;
            }
            maxAbs = std::max(maxAbs, d);
            if(xRef->val[j] != 0.0)
            {
                maxRel = std::max(maxRel, d/fabs(xRef->val[j]));
            }
        }
        printf("BCSR vs point CRS (original numbering): max abs dev = %e, max rel dev = %e, bit-identical %lld of %lld (%.4f %%)\n",
               maxAbs, maxRel, nExact, nTotal, 100.0*nExact/(double)nTotal);
        findMaxDeviations(xRef, xMapped);
        bool flag = checkEqual(xRef, xMapped, param.tol);
        if(!flag)
        {
            printf("Power calculation failed\n");
        }
        else
        {
            printf("Power calculation success\n");
        }
        delete xRef;
        delete xMapped;
        delete pmat;
    }

    delete mat;
    delete xRACE;
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif
    return 0;
}
