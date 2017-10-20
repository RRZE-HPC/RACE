#include "densemat.h"
#include "math.h"

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

bool checkEqual(const densemat* lhs, const densemat* rhs, double tol)
{
    if(lhs->nrows != rhs->nrows)
    {
        printf("Densemat dimension differs\n");
        return false;
    }

    int nrows = lhs->nrows;
    for(int row=0; row<nrows; ++row)
    {
        if( fabs((lhs->val[row]-rhs->val[row])/lhs->val[row]) > tol )
        {
            printf("Densemat deviation @ idx %d lhs = %f, rhs = %f\n", row, lhs->val[row], rhs->val[row]);
            return false;
        }
    }

    return true;
}
