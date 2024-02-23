#include "npy_cblas.h"

#define WRAP_F77(a) a##_

/* This intermediate C-file is necessary also for the Fortran wrapper:
   - The cblas_* functions cannot be called from Fortran directly
     as they do not take all of their arguments as pointers
     (as is done in Fortran).
*/
void WRAP_F77(acc_cdotc_sub)(const int *N, const void *X, const int *incX,
			     const void *Y, const int *incY, void *dotc)
{
    cblas_cdotc_sub(*N, X, *incX, Y, *incY, dotc);
}

void WRAP_F77(acc_cdotu_sub)(const int *N, const void *X, const int *incX,
			     const void *Y, const int *incY, void *dotu)
{
    cblas_cdotu_sub(*N, X, *incX, Y, *incY, dotu);
}

void WRAP_F77(acc_zdotc_sub)(const int *N, const void *X, const int *incX,
			     const void *Y, const int *incY, void *dotc)
{
    cblas_zdotc_sub(*N, X, *incX, Y, *incY, dotc);
}

void WRAP_F77(acc_zdotu_sub)(const int *N, const void *X, const int *incX,
			     const void *Y, const int *incY, void *dotu)
{
    cblas_zdotu_sub(*N, X, *incX, Y, *incY, dotu);
}
