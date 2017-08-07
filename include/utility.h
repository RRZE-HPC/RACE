#ifndef NAME_UTILITY_H
#define NAME_UTILITY_H

#include <algorithm>
#include <iterator>


template <typename T> inline void sort(T *arr, int range_lo, int range_hi, bool rev=false)
{
    if(rev == false) {
        std::stable_sort(arr+range_lo,arr+range_hi);
    } else {
        std::stable_sort(arr+range_lo,arr+range_hi, [&](const int& a, const int& b) { return (arr[a] > arr[b]); });
    }
}


template <typename T> inline void sortPerm(T *arr, int *perm, int range_lo, int range_hi, bool rev=false)
{
    if(rev == false) {
        std::stable_sort(perm+range_lo, perm+range_hi, [&](const int& a, const int& b) {return (arr[a] < arr[b]); });
    } else {
        std::stable_sort(perm+range_lo, perm+range_hi, [&](const int& a, const int& b) {return (arr[a] > arr[b]); });
    }
}

template <typename T> inline T sumArr(T *arr, int len)
{
    T sumVal = 0;
    for(int i=0; i<len; ++i)
    {
        sumVal += arr[i];
    }
    return sumVal;
}

template <typename T> inline T maxArr(T *arr, int len)
{
    T maxVal = 0;
    for(int i=0; i<len; ++i)
    {
        if(arr[i]>maxVal)
        {
            maxVal = arr[i];
        }
    }

    return maxVal;
}


//updates first prmutation array based on the current permutation
inline void updatePerm(int **mainPerm, int *currPerm, int len)
{
    int *totPerm = new int [len];

    if(*mainPerm)
    {
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<len; ++i)
        {
            totPerm[i] = (*mainPerm)[currPerm[i]];
        }

    }
    else
    {
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<len; ++i)
        {
            totPerm[i] = currPerm[i];
        }

    }

    //swap mainPerm and totPerm
    int *temp = (*mainPerm);
    (*mainPerm) = totPerm;
    totPerm = temp;

    delete[] totPerm;
}

void DUMMY(int ctr, bool flag);


#endif
