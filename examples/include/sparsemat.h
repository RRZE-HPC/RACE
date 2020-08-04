#ifndef RACE_SPARSEMAT_H
#define RACE_SPARSEMAT_H

#include <RACE/interface.h>
#include <algorithm>
#include <iterator>

using namespace RACE;

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
    Interface* ce;
    int *rowPtr, *col;
    double *val;

/*    int nrows_bcsr, nnz_bcsr;
    int *rowPtr_bcsr, *col_bcsr;
    double *val_bcsr;*/
    int block_size;

    int *rcmPerm, *rcmInvPerm;

    bool readFile(char* filename);
    bool convertToBCSR(int b_r);
    bool writeFile(char* filename);
    void makeDiagFirst();
    void doRCM();
    void doRCMPermute();
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

    sparsemat();
    ~sparsemat();

};

//Used for NUMA aware allocation, sometimes just
//using first touch policy doesn't give the best performance
struct NUMAmat
{
    sparsemat *mat;
    int *nrows, *nnz, **rowPtr, **col;
    double **val;
    int NUMAdomains;
    std::vector<int> splitRows;

    NUMAmat(sparsemat *mat_, bool manual=false, std::vector<int> splitRows_={-2});
    ~NUMAmat();

    private:
    std::vector<int> getRACEPowerSplit();
};

#endif
