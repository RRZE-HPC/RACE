#ifndef RACE_SPARSEMAT_H
#define RACE_SPARSEMAT_H

#include <RACE/interface.h>
#include <algorithm>
#include <iterator>


template <typename T> void sort_perm(T *arr, int *perm, int len, bool rev=false)
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
    RACEInterface* ce;
    int *rowPtr, *col;
    double *val;

    bool readFile(char* filename);
    bool writeFile(char* filename);
    void makeDiagFirst();
    void colorAndPermute(dist_t dist, int nthreads, int smt=1, PinMethod pinMethod=FILL);
    void permute(int* perm, int* invPerm);
    void NUMA_init(bool symmPart);
    void pinOMP(int nthreads);

    /* For symmetric computations */
    int nnz_symm;
    int *rowPtr_symm, *col_symm;
    double *val_symm;
    bool computeSymmData();

    sparsemat();
    ~sparsemat();

};

#endif
