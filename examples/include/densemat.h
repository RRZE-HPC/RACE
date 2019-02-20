#ifndef RACE_DENSEMAT_H
#define RACE_DENSEMAT_H

#include <functional>

struct densemat
{
    int nrows;
    int ncols;
    bool viewMat;
    double *val;

    void setVal(double value);
    void setRand();
    void setFn(std::function<double(int)> fn);
    void setFn(std::function<double(void)> fn);
    void print();
    densemat* view(int start_col, int end_col);
    densemat(int nrows, int ncols=1, bool view=false);
    ~densemat();

};

bool checkEqual(const densemat* lhs, const densemat* rhs, double tol=1e-4);

#endif
