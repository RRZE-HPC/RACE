#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include "densemat.h"
#include "sparsemat.h"
#include "generator.h"

#define BUILD_CMD(kernelFile, outFile, veclen, nblocks)\
    char* buildCmd;\
    asprintf(&buildCmd, "sed -e 's@#define VECLEN 4@#define VECLEN %d@g' %s > %s && sed -i 's@#define NCOLS 4@#define NCOLS %d@g' %s && make genKernel", veclen, kernelFile, outFile, nblocks, outFile);\
    system(buildCmd);\
    free(buildCmd);\

generator::generator(int veclen, int nblocks)
{
    libHandle = dlopen("libgenKernel.so", RTLD_LOCAL|RTLD_NOW);
    dlclose(libHandle);
    BUILD_CMD("../src/tmpl/kernels.cpp", "../src/gen/kernels.cpp", veclen, nblocks);
    libHandle = dlopen("libgenKernel.so", RTLD_LOCAL|RTLD_NOW);
    if(!libHandle)
    {
        printf("Could not open generated file\n");
    }
}

fn_t generator::getFn(char* name)
{
    dlerror();
    fn_t dynKernel = (fn_t) dlsym(libHandle, name);

    if(dlerror())
    {
        printf("Could not find %s in symbolic table\n", name);
    }

    return dynKernel;
}
