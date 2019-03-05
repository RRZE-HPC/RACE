#include "densemat.h"
#include <cmath>

//currently in column major order
densemat::densemat(int nrows_, int ncols_, bool view_):nrows(nrows_), ncols(ncols_), viewMat(view_)
{
    if(!viewMat)
    {
        val = new double[nrows*ncols];

        for(int j=0; j<ncols; ++j)
        {
#pragma omp parallel for schedule(static)
            for(int i=0; i<nrows; ++i)
            {
                val[j*nrows+i] = 0.0;
            }
        }
    }
}

densemat::~densemat()
{
    if(val && !viewMat)
    {
        delete[] val;
        val = NULL;
    }
}

//view only to columns, can't split rows
densemat* densemat::view(int start_col, int end_col)
{
    densemat* newView = new densemat(this->nrows, (end_col-start_col)+1, true);
    newView->val = &(this->val[start_col*nrows]);
    return newView;
}

void densemat::setVal(double value)
{
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            val[j*nrows+i] = value;
        }
    }
}

void densemat::setFn(std::function<double(int)> fn)
{
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            val[j*nrows+i] = fn(i);
        }
    }
}

void densemat::setFn(std::function<double(void)> fn)
{
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            val[j*nrows+i] = fn();
        }
    }
}

void densemat::setRand()
{
    setFn(rand);
}

bool checkEqual(const densemat* lhs, const densemat* rhs, double tol)
{
    if((lhs->nrows != rhs->nrows) || (lhs->ncols != rhs->ncols))
    {
        printf("Densemat dimension differs\n");
        return false;
    }

    int nrows = lhs->nrows;
    int ncols = lhs->ncols;

    for(int col=0; col<ncols; ++col)
    {
        for(int row=0; row<nrows; ++row)
        {
            double dev = fabs((lhs->val[col*nrows+row]-rhs->val[col*nrows+row]));//lhs->val[col*nrows+row]);
            if( dev > tol )
            {
                printf("Densemat deviation @ idx (%d, %d) lhs = %.18f, rhs = %.18f, dev = %.18f\n", row, col, lhs->val[col*nrows+row], rhs->val[col*nrows+row], dev);
                return false;
            }
        }
    }

    return true;
}

void densemat::print()
{
    for(int row=0; row<nrows; ++row)
    {
        for(int col=0; col<ncols; ++col)
        {
            printf("%f  ", val[col*nrows+row]);
        }
        printf("\n");
    }
}
