#ifndef RACE_DENSEMAT_H
#define RACE_DENSEMAT_H

#include <functional>
#include <vector>

struct densemat
{
    int nrows;
    int ncols;
    bool viewMat;
    double *val;

    void setVal(double value);
    void copyVal(densemat *src);
    void setRand();
    void setFn(std::function<double(int)> fn);
    void setFn(std::function<double(void)> fn);
    void axpby(densemat* x, densemat* y, double a, double b);
    std::vector<double> dot(densemat* x);
    void print();
    densemat* view(int start_col, int end_col);
    densemat(int nrows, int ncols=1, bool view=false);
    ~densemat();
    void initCover(int nrows_, int ncols_, double* val_);
};

bool findMaxDeviations(const densemat* lhs, const densemat* rhs);
bool checkEqual(const densemat* lhs, const densemat* rhs, double tol=1e-4, bool relative=false);

#endif
