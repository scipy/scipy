#ifndef NNLS_H
#define NNLS_H

#include <math.h>
#include "npy_cblas.h"
#include "blaslapack_declarations.h"

void __nnls(const CBLAS_INT m, const CBLAS_INT n, double* restrict a, double* restrict b, double* restrict x, double* restrict w, double* restrict zz, CBLAS_INT* restrict indices, const Py_ssize_t maxiter, double* rnorm, int* info);


#endif
