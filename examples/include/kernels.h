#ifndef RACE_KERNELS_H
#define RACE_KERNELS_H

#include "sparsemat.h"
#include "densemat.h"
#include <complex>
#include <mkl.h>

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

struct kernelArg_split
{
    sparsemat* mat;
    densemat* a;
    densemat* b;
    densemat* x;
};


void spmv(sparsemat* mat, densemat* x);
void spmv(densemat* b, sparsemat* mat, densemat* x);
void plain_spmv(densemat* b, sparsemat* mat, densemat* x);
void spmtv(densemat* b, sparsemat* mat, densemat* x);
void gs(densemat* b, sparsemat* mat, densemat* x);
void gs_serial(densemat* b, sparsemat* mat, densemat* x);
void kacz(densemat* b, sparsemat* mat, densemat* x);
void kacz_serial(densemat* b, sparsemat* mat, densemat* x);
void symm_spmv(densemat* b, sparsemat* mat, densemat* x);
void plain_spmv(sparsemat* mat, densemat* x);
void plain_spmv_only_highest(sparsemat* mat, densemat* x, int power);
void plain_spmv_numa(NUMAmat* mat, densemat* x);


void mkl_spmv(densemat* b, sparsemat* mat, densemat* x, bool symm=false);
void mkl_spmv(sparsemat* mat, densemat* x, bool symm=false);
sparse_matrix_t* mkl_ie_setup(sparsemat* mat, int niter, bool symm=false);
void mkl_ie_spmv(densemat* b, sparse_matrix_t* mat, densemat* x, bool symm=false);
void mkl_ie_spmv(sparse_matrix_t* mat, densemat* x, bool symm=false);
void mkl_ie_free(sparse_matrix_t* A);
void matPower(sparsemat* A, int power, densemat *x);
void matPower_only_highest(sparsemat* A, int power, densemat *x);
void matPowerNuma(NUMAmat* A, int power, densemat *x);
void matPowerBCSR(sparsemat* A, int power, densemat *x);

struct splitHandle
{
    sparsemat* U;
    sparsemat* L;
    densemat* Un_x; //tmp workspace
    densemat* Ln_x; //tmp workspace

    splitHandle();
};

splitHandle* matPower_split_init(sparsemat* L, sparsemat *U);
void matPower_split_destroy(splitHandle* handle);
void matPower_split(splitHandle *handle, int power, densemat *x);

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

//convenience macros
#define ENCODE_TO_VOID_SPLIT(mat_en, a_en, b_en, x_en, _NAME_)\
    kernelArg_split *arg_encode_##_NAME_ = new kernelArg_split;\
    arg_encode_##_NAME_->mat = mat_en;\
    arg_encode_##_NAME_->a = a_en;\
    arg_encode_##_NAME_->b = b_en;\
    arg_encode_##_NAME_->x = x_en;\
    void* voidArg_##_NAME_ = (void*) arg_encode_##_NAME_;\

#define UPDATE_X_SPLIT(x_en, _NAME_)\
    arg_encode_##_NAME_->x = x_en;\

#define DELETE_ARG_SPLIT(_NAME_)\
    delete arg_encode_##_NAME_;\


#define DECODE_FROM_VOID_SPLIT(voidArg, _NAME_)\
    kernelArg_split* arg_decode_##_NAME_ = (kernelArg_split*) voidArg;\
    sparsemat* mat = arg_decode_##_NAME_->mat;\
    densemat* a = arg_decode_##_NAME_->a;\
    densemat* b = arg_decode_##_NAME_->b;\
    densemat* x = arg_decode_##_NAME_->x;\


#endif
