#ifndef __NNLS_H
#define __NNLS_H
#include <math.h>

double ddot_(int* n, double* dx, int* incx, double* dy, int* incy);
void dlarf_(char* side, int* m, int* n, double* v, int* incv, double* tau, double* c, int* ldc, double* work);
void dlarfgp_(int* n, double* alpha, double* x, int* incx, double* tau);
void dlartgp_(double* f, double* g, double* cs, double* sn, double* r);
double dnrm2_(int* n, double* x, int* incx);

void
__nnls(const int m, const int n, double* restrict a, double* restrict b,
       double* restrict x, double* restrict w, double* restrict zz,
       int* restrict indices, const int maxiter, double* rnorm, int* info);


#endif
