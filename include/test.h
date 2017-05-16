#ifndef NAME_TEST_H
#define NAME_TEST_H

inline void test(double *a, double *b, double *c, double *d, int start, int end, int ctr)
{
/*    int len = static_cast<int>(10*1024*1024/8.0);
    double *a = new double[len];
    double *b = new double[len];
    double *c = new double[len];
    double *d = new double[len];
*/
    for(int time=0; time<ctr; ++time)
    {
        for(int i=start; i<end; ++i)
        {
            a[i] = b[i]+c[i]*d[i];
        }
    }
/*
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
*/
}
#endif
