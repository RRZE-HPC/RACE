#include "mmio.h"
#include "time.h"
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <unistd.h>

#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif
#include "densemat.h"
#include "kernels.h"
#include "parse.h"
#include "sparsemat.h"
#include "timer.h"
#include <algorithm>
#include <iostream>

// #define VALIDATE_wo_PERM

// Number of times the whole benchmark loop is repeated; the *best* (minimum)
// synchronized time is reported. 3 is usually enough to kill most of the
// run-to-run jitter you see in the heatmaps. Override with -DNUM_REPEATS=n.
#ifndef NUM_REPEATS
#define NUM_REPEATS 3
#endif

// If defined, an MPI_Barrier is issued after every single power iteration.
// This removes drift between ranks at the cost of adding barrier overhead to
// the measurement -- useful for debugging imbalance, not for peak numbers.
// #define SYNC_EACH_ITER

static int myRank = 0;
static int nRanks = 1;

#define RANK0_PRINT(...)                                                       \
  do {                                                                         \
    if (myRank == 0) {                                                         \
      printf(__VA_ARGS__);                                                     \
      fflush(stdout);                                                          \
    }                                                                          \
  } while (0)

struct timeStat {
  double max; // slowest rank -> this is the "fair" time
  double min; // fastest rank
  double avg;
  double imbalance; // (max-min)/max in percent
};

// Collective: every rank contributes its local time, everybody gets the same
// statistics back so that any rank can compute the same performance number.
static timeStat reduceTime(double localTime) {
  timeStat s;
  double sum = 0.0;

  MPI_Allreduce(&localTime, &s.max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&localTime, &s.min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&localTime, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  s.avg = sum / (double)nRanks;
  s.imbalance = (s.max > 0.0) ? 100.0 * (s.max - s.min) / s.max : 0.0;

  return s;
}

// flops here is the *per rank* GFlop count.
static void reportPerf(const char *name, timeStat t, double flops) {
  RANK0_PRINT("%-16s : time(max) = %8.5f s | time(min) = %8.5f s | time(avg) = "
              "%8.5f s | imbalance = %5.2f %%\n",
              name, t.max, t.min, t.avg, t.imbalance);
  RANK0_PRINT("%-16s : perf/rank = %8.4f GFlop/s | aggregate = %8.4f GFlop/s "
              "(%d ranks)\n",
              name, flops / t.max, nRanks * flops / t.max, nRanks);
}

void capitalize(char *beg) {
  int i = 0;
  while (beg[i]) {
    beg[i] = toupper(beg[i]);
    ++i;
  }
}

#define PERF_RUN(kernel, flopPerNnz)                                           \
  {                                                                            \
    int iter = param.iter;                                                     \
    double time = 0;                                                           \
    double nnz_update = ((double)mat->nnz) * iter * 1e-9;                      \
    sleep(1);                                                                  \
    MPI_Barrier(MPI_COMM_WORLD);                                               \
    INIT_TIMER(kernel);                                                        \
    START_TIMER(kernel);                                                       \
    kernel(b, mat, x, iter);                                                   \
    STOP_TIMER(kernel);                                                        \
    time = GET_TIMER(kernel);                                                  \
    timeStat ts = reduceTime(time);                                            \
    char *capsKernel;                                                          \
    asprintf(&capsKernel, "%s", #kernel);                                      \
    capitalize(capsKernel);                                                    \
    RANK0_PRINT("%10s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel,         \
                flopPerNnz * nnz_update / (ts.max), ts.max);                   \
    free(capsKernel);                                                          \
  }

int findMetricId(int group_id, std::string toFind) {
#ifdef LIKWID_PERFMON
  int numMetrics = perfmon_getNumberOfMetrics(group_id);

  int dataVol_metric_id = numMetrics - 1;
  for (int i = 0; i < numMetrics; ++i) {
    std::string currMetric(perfmon_getMetricName(group_id, i));
    if (myRank == 0) {
      std::cout << "currMetric = " << currMetric << std::endl;
    }
    if (currMetric.find(toFind) != std::string::npos) {
      dataVol_metric_id = i;
    }
  }

  return dataVol_metric_id;
#else
  return -1;
#endif
}

int main(int argc, char *argv[]) {
  int provided = 0;
  // OpenMP threads inside, only the master thread does MPI -> FUNNELED
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  if (provided < MPI_THREAD_FUNNELED) {
    RANK0_PRINT("WARNING: MPI library provides thread level %d, "
                "MPI_THREAD_FUNNELED (%d) requested\n",
                provided, MPI_THREAD_FUNNELED);
  }

#ifdef LIKWID_PERFMON
  LIKWID_MARKER_INIT;
#endif

  int nthreads = 1;
#pragma omp parallel
  {
#pragma omp single
    nthreads = omp_get_num_threads();
  }
  RANK0_PRINT("Running with %d MPI rank(s) x %d OpenMP thread(s)\n", nRanks,
              nthreads);

  parser param;
  if (!param.parse_arg(argc, argv)) {
    RANK0_PRINT("Error in reading parameters\n");
  }

  sparsemat *mat = new sparsemat;

  RANK0_PRINT("Reading matrix file\n");
  // Stagger the I/O a bit so that N ranks don't hammer the same file at once.
  for (int r = 0; r < nRanks; ++r) {
    if (r == myRank) {
      if (!mat->readFile(param.mat_file)) {
        printf("[rank %d] Error in reading sparse matrix file\n", myRank);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  int NROWS = mat->nrows;
  bool randInit = false;
  double initVal = 1 / (double)NROWS;
  densemat *xRAND = NULL;
  if (randInit) {
    xRAND = new densemat(NROWS);
    xRAND->setRand();
  }
  int power = param.iter;
  RANK0_PRINT("power = %d\n", power);

  densemat *xTRAD = NULL;
#ifdef VALIDATE_wo_PERM
  if (param.validate) {
    xTRAD = new densemat(NROWS, power + 1);
    densemat *xTRAD_0 = xTRAD->view(0, 0);

    if (randInit) {
      xTRAD_0->copyVal(xRAND);
    } else {
      xTRAD_0->setVal(initVal);
    }

    for (int pow = 0; pow < power; ++pow) {
      densemat *x = xTRAD->view(pow, pow);
      plain_spmv(mat, x);
    }
  }
#endif

  RANK0_PRINT("Preparing matrix for power calculation\n");

  MPI_Barrier(MPI_COMM_WORLD);
  INIT_TIMER(pre_process);
  START_TIMER(pre_process);
  if (param.RCM_flag) {
    mat->doRCM();
  }
  mat->prepareForPower(power, param.cache_size, param.cores, param.smt,
                       param.pin);
  STOP_TIMER(pre_process);
  timeStat preStat = reduceTime(GET_TIMER(pre_process));
  RANK0_PRINT("Total pre-processing time (max over ranks) = %f s\n",
              preStat.max);

  INFO_PRINT("Matrix statistics");
  INFO_PRINT("Nrows = %d, NNZ = %d, NNZR = %f\n", mat->nrows, mat->nnz,
             mat->nnz / ((double)mat->nrows));

  densemat *xRACE = new densemat(NROWS, power + 1);
  densemat *xRACE_0 = xRACE->view(0, 0);
  if (randInit) {
    xRACE_0->copyVal(xRAND);
  } else {
    xRACE_0->setVal(initVal);
  }

  RANK0_PRINT("calculation started\n");

  //-------------------------------------------------------------------------
  // Determine the number of iterations. This MUST be identical on all ranks,
  // otherwise the barriers below compare apples with oranges. We therefore
  // base it on the slowest rank and broadcast it.
  //-------------------------------------------------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  INIT_TIMER(matPower_init);
  START_TIMER(matPower_init);
  for (int iter = 0; iter < 10; ++iter) {
    matPower(mat, power, xRACE);
  }
  STOP_TIMER(matPower_init);
  timeStat initStat = reduceTime(GET_TIMER(matPower_init));

  int iterations = std::max(1, (int)(1.2 * 10 / initStat.max));
  MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
  RANK0_PRINT("Num iterations =  %d (calibration time = %f s)\n", iterations,
              initStat.max);

  // GFlop performed by ONE rank
  double flops = 2.0 * power * iterations * (double)mat->nnz * 1e-9;

  //-------------------------------------------------------------------------
  // Reference: plain SpMV power
  //-------------------------------------------------------------------------
  double spmvBestTime = 1e300;
  timeStat spmvStat = {0, 0, 0, 0};
  if (param.validate) {
    densemat *xTRAD_perf = new densemat(NROWS, power + 1);
#ifndef VALIDATE_wo_PERM
    xTRAD = xTRAD_perf;
#endif
    densemat *xTRAD_0 = xTRAD_perf->view(0, 0);

    for (int rep = 0; rep < NUM_REPEATS; ++rep) {
      xTRAD_perf->setVal(0);
      if (randInit) {
        xTRAD_0->copyVal(xRAND);
      } else {
        xTRAD_0->setVal(initVal);
      }

      sleep(1);
      MPI_Barrier(MPI_COMM_WORLD); // everybody starts together
      INIT_TIMER(spmvPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
      {
        LIKWID_MARKER_START("SpMV_power");
      }
#endif
      START_TIMER(spmvPower);
      for (int iter = 0; iter < iterations; ++iter) {
        for (int pow = 0; pow < power; ++pow) {
          densemat *x = xTRAD_perf->view(pow, pow);
          plain_spmv(mat, x);
        }
#ifdef SYNC_EACH_ITER
        MPI_Barrier(MPI_COMM_WORLD);
#endif
      }
      STOP_TIMER(spmvPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
      {
        LIKWID_MARKER_STOP("SpMV_power");
      }
#endif
      timeStat s = reduceTime(GET_TIMER(spmvPower));
      RANK0_PRINT("  [SpMV  rep %d/%d] max = %8.5f s, imbalance = %5.2f %%\n",
                  rep + 1, NUM_REPEATS, s.max, s.imbalance);
      if (s.max < spmvBestTime) {
        spmvBestTime = s.max;
        spmvStat = s;
      }
    }

    reportPerf("SpMV_power", spmvStat, flops);

#ifdef VALIDATE_wo_PERM
    delete xTRAD_perf;
#endif
    sleep(1);
  }

  //-------------------------------------------------------------------------
  // RACE matrix power kernel
  //-------------------------------------------------------------------------
  double raceBestTime = 1e300;
  timeStat raceStat = {0, 0, 0, 0};

  for (int rep = 0; rep < NUM_REPEATS; ++rep) {
    xRACE->setVal(0);
    if (randInit) {
      xRACE_0->copyVal(xRAND);
    } else {
      xRACE_0->setVal(initVal);
    }

    sleep(1);
    MPI_Barrier(MPI_COMM_WORLD); // fair start for all ranks
    INIT_TIMER(matPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
    {
      LIKWID_MARKER_START("RACE_power");
    }
#endif
    START_TIMER(matPower);
    for (int iter = 0; iter < iterations; ++iter) {
      matPower(mat, power, xRACE);
#ifdef SYNC_EACH_ITER
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    STOP_TIMER(matPower);
#ifdef LIKWID_PERFMON
#pragma omp parallel
    {
      LIKWID_MARKER_STOP("RACE_power");
    }
#endif
    timeStat s = reduceTime(GET_TIMER(matPower));
    RANK0_PRINT("  [RACE  rep %d/%d] max = %8.5f s, imbalance = %5.2f %%\n",
                rep + 1, NUM_REPEATS, s.max, s.imbalance);
    if (s.max < raceBestTime) {
      raceBestTime = s.max;
      raceStat = s;
    }
  }

  reportPerf("RACE_power", raceStat, flops);

  if (param.validate && spmvBestTime < 1e299) {
    RANK0_PRINT("Speedup (RACE/SpMV, synchronized best times) = %6.3f\n",
                spmvBestTime / raceBestTime);
  }

  //-------------------------------------------------------------------------
  // Validation
  //-------------------------------------------------------------------------
  if (param.validate) {
#ifdef VALIDATE_wo_PERM
    densemat *xTRAD_permuted = mat->permute_densemat(xTRAD);
#else
    densemat *xTRAD_permuted = xTRAD;
#endif

    xRACE->setVal(0);
    if (randInit) {
      xRACE_0->copyVal(xRAND);
    } else {
      xRACE_0->setVal(initVal);
    }

#ifdef VALIDATE_wo_PERM
    densemat *xRACE_permuted = mat->permute_densemat(xRACE);
#else
    densemat *xRACE_permuted = xRACE;
#endif

    matPower(mat, power, xRACE_permuted);

    findMaxDeviations(xTRAD_permuted, xRACE_permuted);
    int localOk = checkEqual(xTRAD_permuted, xRACE_permuted, param.tol) ? 1 : 0;

    int globalOk = 0;
    MPI_Allreduce(&localOk, &globalOk, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    if (!localOk) {
      printf("[rank %d] Power calculation failed\n", myRank);
    }
    RANK0_PRINT("Power calculation %s (all ranks)\n",
                globalOk ? "success" : "FAILED");

#ifdef VALIDATE_wo_PERM
    delete xTRAD_permuted;
    delete xRACE_permuted;
#endif
    delete xTRAD;
  }

  delete mat;
  delete xRACE;
  if (randInit) {
    delete xRAND;
  }

#ifdef LIKWID_PERFMON
  LIKWID_MARKER_CLOSE;
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}