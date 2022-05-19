/*
 * =======================================================================================
 *
 *   RACE: Recursicve Algebraic Coloring Engine
 *   Copyright (C) 2019, RRZE, Friedrich-Alexander-Universität Erlangen-Nürnberg
 *   Author: Christie Alappat
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * =======================================================================================
 */

#ifndef RACE_UTILITY_H
#define RACE_UTILITY_H

#include <algorithm>
#include <iterator>
#include <vector>
#include <unistd.h>
#include <sys/mman.h>
#include <stdio.h>
#include <string.h>
#include <cmath>

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
        char *copyEnvVal = (char*) malloc(strlen(envVal) + 1);
        strcpy(copyEnvVal, envVal);

        char* token = strtok(copyEnvVal, ",");
        while(token != NULL)
        {
            val_vec.push_back(atoi(token));
            token = strtok(NULL, ",");
        }

        free(copyEnvVal);
    }
}

//updates first permutation array based on the current permutation
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


//updates first permutation array based on the current permutation
inline void updatePerm(int **mainPerm, std::vector<int> currPerm, int len, int fullLen)
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

inline int workingBoundaryLength_base(int highestPower)
{
    return static_cast<int>(ceil(highestPower/2.0))-1;
}

#endif
