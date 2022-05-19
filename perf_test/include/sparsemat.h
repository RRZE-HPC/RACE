#ifndef RACE_SPARSEMAT_H
#define RACE_SPARSEMAT_H

#include "interface.h"
#include <algorithm>
#include <iterator>

using namespace RACE;

template <typename T> inline void sort_perm(T *arr, int *perm, int len, bool rev=false)
{
    if(rev == false) {
        std::stable_sort(perm+0, perm+len, [&](const int& a, const int& b) {return (arr[a] < arr[b]); });
    } else {
        std::stable_sort(perm+0, perm+len, [&](const int& a, const int& b) {return (arr[a] > arr[b]); });
    }
}


struct sparsemat
{
    int nrows, nnz;
    //interface to coloring engine
    Interface* ce;
    double *val;
    int *rowPtr, *col;

    bool readFile(char* filename);
    bool writeFile(char* filename);
    void makeDiagFirst();
    void doRCM();
    int prepareForPower(int highestPower, int numSharedCache, double cacheSize, int nthreads, int smt=1, PinMethod pinMethod=FILL);
    int colorAndPermute(dist d, int nthreads, int smt=1, PinMethod pinMethod=FILL);
    double colorEff();
    int maxStageDepth();
    void permute(int* perm, int* invPerm, bool RACEalloc=false);
    void numaInit(bool RACEalloc=false);
    void pinOMP(int nthreads);

    /* For symmetric computations */
    int nnz_symm;
    int *rowPtr_symm, *col_symm;
    double *val_symm;
    bool computeSymmData();

    int *rcmPerm, *rcmInvPerm;
    sparsemat();
    ~sparsemat();

};

#endif
