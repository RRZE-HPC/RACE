#ifndef RACE_UTILITY_H
#define RACE_UTILITY_H

#include <algorithm>
#include <iterator>
#include <vector>
#include <unistd.h>
#include <sys/mman.h>

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

template <typename T> inline void getEnv(std::string envName, std::vector<T>& val_vec)
{
    char *envVal = getenv(envName.c_str());

    if(envVal != NULL)
    {
        char* token = strtok(envVal, ",");
        while(token != NULL)
        {
            val_vec.push_back(atoi(token));
            token = strtok(NULL, ",");
        }
    }
}

//updates first peermutation array based on the current permutation
inline void updatePerm(int **mainPerm, int *currPerm, int len, int fullLen)
{
    int *totPerm = new int [fullLen];

    if(*mainPerm)
    {
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<len; ++i)
        {
            totPerm[i] = (*mainPerm)[currPerm[i]];
        }
        for(int i=len; i<fullLen; ++i)
        {
            totPerm[i] = (*mainPerm)[i];
        }
    }
    else
    {
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<len; ++i)
        {
            totPerm[i] = currPerm[i];
        }
        for(int i=len; i<fullLen; ++i)
        {
            totPerm[i] = i;
        }
    }

    //swap mainPerm and totPerm
    int *temp = (*mainPerm);
    (*mainPerm) = totPerm;
    totPerm = temp;

    delete[] totPerm;
}

void DUMMY(int ctr, bool flag);

void DUMMY_spin(volatile unsigned int *a);

inline int roundDouble(double val)
{
    int val_i = static_cast<int>(val);
    return ( (((val - val_i))<0.5)?val_i:(val_i+1) );
}

//Taken from :
//https://www3.physnet.uni-hamburg.de/physnet/Tru64-Unix/HTML/APS33DTE/DOCU_005.HTM
inline int lock_memory(char   *addr,
        size_t  size)
{
    unsigned long    page_offset, page_size;

    page_size = sysconf(_SC_PAGE_SIZE);
    page_offset = (unsigned long) addr % page_size;

    addr -= page_offset;  /* Adjust addr to page boundary */
    size += page_offset;  /* Adjust size with page_offset */

    return ( mlock(addr, size) );  /* Lock the memory */
}

inline int unlock_memory(char   *addr,
        size_t  size)
{
    unsigned long    page_offset, page_size;

    page_size = sysconf(_SC_PAGE_SIZE);
    page_offset = (unsigned long) addr % page_size;

    addr -= page_offset;  /* Adjust addr to page boundary */
    size += page_offset;  /* Adjust size with page_offset */

    return ( munlock(addr, size) );  /* Unlock the memory */
}

#endif
