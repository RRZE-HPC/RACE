#include "densemat.h"

densemat::densemat(int nrows_):nrows(nrows_)
{
    val = new double[nrows];

#pragma omp parallel for
    for(int i=0; i<nrows; ++i)
    {
        val[i] = 0.0;
    }
}

densemat::~densemat()
{
    if(val)
    {
        delete[] val;
    }
}

void densemat::setVal(double value)
{
#pragma omp parallel for
    for(int i=0; i<nrows; ++i)
    {
        val[i] = value;
    }
}

void densemat::setFn(std::function<double(int)> fn)
{
#pragma omp parallel for
    for(int i=0; i<nrows; ++i)
    {
        val[i] = fn(i);
    }
}

void densemat::setFn(std::function<double(void)> fn)
{
#pragma omp parallel for
    for(int i=0; i<nrows; ++i)
    {
        val[i] = fn();
    }
}

void densemat::setRand()
{
    setFn(rand);
}

