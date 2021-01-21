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
#include <map>
#include <iostream>

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

    int NROWS = mat->nrows;

    std::map<int, int> nnzr_bucket;
    for(int row=0; row<NROWS; ++row)
    {
        int curr_nnzr = (mat->rowPtr[row+1]-mat->rowPtr[row]);
        nnzr_bucket[curr_nnzr] += 1;
    }

    printf("nnzr -> freq\n");
    for(auto it = nnzr_bucket.cbegin(); it != nnzr_bucket.cend(); ++it)
    {
        std::cout<<it->first<<"->"<<it->second<<std::endl;
    }

    printf("Nrows = %d\n", NROWS);
    printf("NNZ = %d\n", mat->nnz);
    delete mat;
}
