#ifndef RACE_DENSEMAT_H
#define RACE_DENSEMAT_H

#include <functional>

struct densemat
{
    int nrows;
    int ncols;
    double *val;
    bool complex_value;
    void setVal(double value, double value_im=0);
    void setRand();
    void setFn(std::function<double(int)> fn);
    void setFn(std::function<double(void)> fn);

    densemat(int nrows, bool complex_value_=false, int ncols=1);
    ~densemat();

};

bool checkEqual(const densemat* lhs, const densemat* rhs, double tol=1e-4, bool complex_value=false);

#endif
