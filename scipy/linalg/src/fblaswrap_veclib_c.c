#include <vecLib/vecLib.h>

//#define WRAP_F77(a) wcblas_##a##_
#define WRAP_F77(a) w##a##_
void WRAP_F77(cdotc)(complex *dotc, const int *N, const complex *X, const int *incX, const complex *Y, const int *incY)
{
    cblas_cdotc_sub(*N, X, *incX, Y, *incY, dotc);
}

void WRAP_F77(cdotu)(complex* dotu, const int *N, const complex *X, const int *incX, const complex *Y, const int *incY)
{
    cblas_cdotu_sub(*N, X, *incX, Y, *incY, dotu);
}

void WRAP_F77(zdotc)(double complex *dotu, const int *N, const double complex *X, const int *incX, const double complex *Y, const int *incY)
{
    cblas_zdotc_sub(*N, X, *incX, Y, *incY, dotu);
}
void WRAP_F77(zdotu)(double complex *dotu, const int *N, const double complex *X, const int *incX, const double complex *Y, const int *incY)
{
    cblas_zdotu_sub(*N, X, *incX, Y, *incY, dotu);
}
