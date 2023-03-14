#include "densemat.h"
#include <cmath>
#include <stdio.h>
#include <omp.h>

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

//similar to view but can pass value too
void densemat::initCover(int nrows_, int ncols_, double* val_)
{
    nrows = nrows_;
    ncols = ncols_;
    val = val_;
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

void densemat::copyVal(densemat* src)
{
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel for
        for(int i=0; i<nrows; ++i)
        {
            val[j*nrows+i] = src->val[j*nrows+i];
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

//this = a*x[] + b*y[]
void densemat::axpby(densemat* x, densemat* y, double a, double b)
{
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel for schedule(static)
        for(int i=0; i<nrows; ++i)
        {
            val[j*nrows+i] = a*x->val[j*nrows+i] + b*y->val[j*nrows+i];
        }
    }
}

//res = dot(this, x)
std::vector<double> densemat::dot(densemat* x)
{
    std::vector<double> dots(ncols);
    for(int j=0; j<ncols; ++j)
    {
        double curDot = 0;
#pragma omp parallel for schedule(static) reduction(+:curDot)
        for(int i=0; i<nrows; ++i)
        {
            curDot += val[j*nrows+i]*x->val[j*nrows+i];
        }
        dots[j] = curDot;
    }
    return dots;
}

void densemat::setRand()
{
    //setFn(rand);
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel
        {
            unsigned int seed = omp_get_thread_num()*100;
#pragma omp for
            for(int i=0; i<nrows; ++i)
            {
                val[j*nrows+i] = rand_r(&seed);
            }
        }
    }

}

//this stops once a deviation higher than tol found
bool checkEqual(const densemat* lhs, const densemat* rhs, double tol, bool relative)
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
            double denom = 1;
            if(relative)
            {
                denom=lhs->val[col*nrows+row];
            }
            double dev = fabs((lhs->val[col*nrows+row]-rhs->val[col*nrows+row])/denom);
            if( dev > tol )
            {
                char* errType = "absolute";
                if(relative)
                {
                    errType = "relative";
                }
                printf("Densemat %s deviation @ idx (%d, %d) lhs = %.18f, rhs = %.18f, dev = %.18f\n", errType, row, col, lhs->val[col*nrows+row], rhs->val[col*nrows+row], dev);
                return false;
            }
        }
    }
    return true;
}


bool findMaxDeviations(const densemat* lhs, const densemat* rhs)
{
    if((lhs->nrows != rhs->nrows) || (lhs->ncols != rhs->ncols))
    {
        printf("Densemat dimension differs\n");
        return false;
    }

    int nrows = lhs->nrows;
    int ncols = lhs->ncols;

    double max_dev = 0;
    double max_rel_dev = 0;

    int row_max_dev=0, col_max_dev=0;
    int row_max_rel_dev=0, col_max_rel_dev=0;

    for(int col=0; col<ncols; ++col)
    {
        for(int row=0; row<nrows; ++row)
        {
            double denom=lhs->val[col*nrows+row];
            double dev = fabs((lhs->val[col*nrows+row]-rhs->val[col*nrows+row]));
            double rel_dev = fabs((lhs->val[col*nrows+row]-rhs->val[col*nrows+row])/denom);
            if(dev > max_dev)
            {
                max_dev = dev;
                row_max_dev = row;
                col_max_dev = col;
            }
            if(rel_dev > max_rel_dev)
            {
                max_rel_dev = rel_dev;
                row_max_rel_dev = row;
                col_max_rel_dev = col;
            }
        }
    }
    printf("Max. densemat deviation @ idx (%d, %d) lhs = %.18f, rhs = %.18f, dev = %.18f\n", row_max_dev, col_max_dev, lhs->val[col_max_dev*nrows+row_max_dev], rhs->val[col_max_dev*nrows+row_max_dev], max_dev);
    printf("Max. relative densemat deviation @ idx (%d, %d) lhs = %.18f, rhs = %.18f, relative dev = %.18f\n", row_max_rel_dev, col_max_rel_dev, lhs->val[col_max_rel_dev*nrows+row_max_rel_dev], rhs->val[col_max_rel_dev*nrows+row_max_rel_dev], max_rel_dev);
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
