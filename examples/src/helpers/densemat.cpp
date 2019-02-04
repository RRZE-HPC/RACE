#include "densemat.h"
#include "math.h"

densemat::densemat(int nrows_, bool complex_value_, int ncols_):nrows(nrows_), complex_value(complex_value_), ncols(ncols_)
{
    int multiple = 1;
    if(!complex_value)
    {
        val = new double[nrows*ncols];
    }
    else
    {
        val = new double[2*nrows*ncols];
        multiple = 2;
    }

#pragma omp parallel for
    for(int i=0; i<multiple*nrows*ncols; ++i)
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

void densemat::setVal(double value, double value_im)
{
    if(!complex_value)
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            for(int j=0; j<ncols; ++j)
            {
                val[i*ncols+j] = value;
            }
        }
    }
    else
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            for(int j=0; j<ncols; ++j)
            {
                val[2*(i*ncols+j)] = value;
                val[2*(i*ncols+j)+1] = value_im;
            }
        }
    }
}

void densemat::setFn(std::function<double(int)> fn)
{
    if(!complex_value)
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            for(int j=0; j<ncols; ++j)
            {
                val[i*ncols+j] = fn(i);
            }
        }
    }
    else
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            for(int j=0; j<ncols; ++j)
            {
                val[2*(i*ncols+j)] = fn(i);
                val[2*(i*ncols+j)+1] = 0;
            }
        }
    }
}

void densemat::setFn(std::function<double(void)> fn)
{
    if(!complex_value)
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            for(int j=0; j<ncols; ++j)
            {
                val[i*ncols+j] = fn();
            }
        }
    }
    else
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            for(int j=0; j<ncols; ++j)
            {
                val[2*(i*ncols+j)] = fn();
                val[2*(i*ncols+j)+1] = 0;
            }
        }
    }
}

void densemat::setRand()
{
    setFn(rand);
}

bool checkEqual(const densemat* lhs, const densemat* rhs, double tol, bool complex_value)
{
    if((lhs->nrows != rhs->nrows) || (lhs->ncols != rhs->ncols))
    {
        printf("Densemat dimension differs\n");
        return false;
    }

    int nrows = lhs->nrows;
    int ncols = lhs->ncols;

    if(!complex_value)
    {
        for(int row=0; row<nrows; ++row)
        {
            for(int j=0; j<ncols; ++j)
            {
                if( fabs((lhs->val[row*ncols+j]-rhs->val[row*ncols+j])/lhs->val[row*ncols+j]) > tol )
                {
                    printf("Densemat deviation @ idx (%d, %d) lhs = %f, rhs = %f\n", row, j, lhs->val[row*ncols+j], rhs->val[row*ncols+j]);
                    return false;
                }
            }
        }
    }
    else
    {
        for(int row=0; row<nrows; ++row)
        {
            for(int j=0; j<ncols; ++j)
            {
                if( fabs( (lhs->val[2*(row*ncols+j)]-rhs->val[2*(row*ncols+j)])/(lhs->val[2*(row*ncols+j)]) ) > (tol) || fabs( (lhs->val[2*(row*ncols+j)+1]-rhs->val[2*(row*ncols+j)+1])/(lhs->val[2*(row*ncols+j)+1]) ) > (tol) )
                {
                    printf("Densemat deviation @ idx (%d, %d) lhs = %f + i%f, rhs = %f + i%f\n", row, j, lhs->val[2*(row*ncols+j)], lhs->val[2*(row*ncols+j)+1], rhs->val[2*(row*ncols+j)], rhs->val[2*(row*ncols+j)+1]);
                    return false;
                }
            }
        }
    }


    return true;
}
