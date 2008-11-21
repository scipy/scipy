#include <vecLib/vecLib.h>

#define WRAP_F77(a) a##_
void WRAP_F77(veclib_cdotc)(const int *N, const complex float *X, const int *incX,
const complex float *Y, const int *incY, complex float *dotc)
{
    cblas_cdotc_sub(*N, X, *incX, Y, *incY, dotc);
}

void WRAP_F77(veclib_cdotu)(const int *N, const complex float *X, const int *incX,
const complex float *Y, const int *incY, complex float* dotu)
{
    cblas_cdotu_sub(*N, X, *incX, Y, *incY, dotu);
}

void WRAP_F77(veclib_zdotc)(const int *N, const double complex *X, const int
*incX, const double complex *Y, const int *incY, double complex *dotu)
{
    cblas_zdotc_sub(*N, X, *incX, Y, *incY, dotu);
}
void WRAP_F77(veclib_zdotu)(const int *N, const double complex *X, const int
*incX, const double complex *Y, const int *incY, double complex *dotu)
{
    cblas_zdotu_sub(*N, X, *incX, Y, *incY, dotu);
}
