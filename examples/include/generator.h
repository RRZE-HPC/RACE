#ifndef _GEN_KERNEL_
#define _GEN_KENREL_

typedef void (*fn_t) (densemat*, sparsemat*, densemat*, int);

fn_t generateKernel(int veclen, int nblocks, char* kernelName);

class generator
{
    public:
        void* libHandle;
        generator(int veclen, int nblocks);
        fn_t getFn(char* name);
};

#endif
