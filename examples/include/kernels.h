#ifndef RACE_KERNELS_H
#define RACE_KERNELS_H

#include "sparsemat.h"
#include "densemat.h"

struct kernelArg
{
    sparsemat* mat;
    densemat* b;
    densemat* x;
};


struct kernelNumaArg
{
    NUMAmat* mat;
    densemat* b;
    densemat* x;
};


void spmv(sparsemat* mat, densemat* x);
void spmv(densemat* b, sparsemat* mat, densemat* x, int iter);
void plain_spmv(densemat* b, sparsemat* mat, densemat* x, int iter);
void spmtv(densemat* b, sparsemat* mat, densemat* x, int iter);
void gs(densemat* b, sparsemat* mat, densemat* x, int iter);
void kacz(densemat* b, sparsemat* mat, densemat* x, int iter);
void symm_spmv(densemat* b, sparsemat* mat, densemat* x, int iter);
void plain_spmv(sparsemat* mat, densemat* x);
void plain_spmv_numa(NUMAmat* mat, densemat* x);

void mkl_spmv(sparsemat* mat, densemat* x);
void matPower(sparsemat* A, int power, densemat *x);
void matPowerNuma(NUMAmat* A, int power, densemat *x);
void matPowerBCSR(sparsemat* A, int power, densemat *x);

//convenience macros
#define ENCODE_TO_VOID(mat_en, b_en, x_en)\
    kernelArg *arg_encode = new kernelArg;\
    arg_encode->mat = mat_en;\
    arg_encode->b = b_en;\
    arg_encode->x = x_en;\
    void* voidArg = (void*) arg_encode;\

#define DELETE_ARG()\
    delete arg_encode;\


#define DECODE_FROM_VOID(voidArg)\
    kernelArg* arg_decode = (kernelArg*) voidArg;\
    sparsemat* mat = arg_decode->mat;\
    densemat* b = arg_decode->b;\
    densemat* x = arg_decode->x;\


#define ENCODE_TO_VOID_NUMA(mat_en, b_en, x_en)\
    kernelNumaArg *arg_encode = new kernelNumaArg;\
    arg_encode->mat = mat_en;\
    arg_encode->b = b_en;\
    arg_encode->x = x_en;\
    void* voidArg = (void*) arg_encode;\


#define DECODE_FROM_VOID_NUMA(voidArg)\
    kernelNumaArg* arg_decode = (kernelNumaArg*) voidArg;\
    NUMAmat* mat = arg_decode->mat;\
    densemat* b = arg_decode->b;\
    densemat* x = arg_decode->x;\


#endif
