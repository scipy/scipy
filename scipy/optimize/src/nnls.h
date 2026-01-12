#ifndef NNLS_H
#define NNLS_H

#include <math.h>
#include "blaslapack_declarations.h"

void __nnls(const int m, const int n, double* restrict a, double* restrict b, double* restrict x, double* restrict w, double* restrict zz, int* restrict indices, const int maxiter, double* rnorm, int* info);


#endif
