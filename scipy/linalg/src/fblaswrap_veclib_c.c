#include <Accelerate/Accelerate.h>

//#define WRAP_F77(a) wcblas_##a##_
#define WRAP_F77(a) w##a##_
void WRAP_F77(cdotc)(void *dotc, const int *N, const void *X, const int *incX, 
		     const void *Y, const int *incY)
{
    cblas_cdotc_sub(*N, X, *incX, Y, *incY, dotc);
}

void WRAP_F77(cdotu)(void *dotu, const int *N, const void *X, const int *incX,
		     const void *Y, const int *incY)
{
    cblas_cdotu_sub(*N, X, *incX, Y, *incY, dotu);
}

void WRAP_F77(zdotc)(void *dotu, const int *N, const void *X, const int *incX,
		     const void *Y, const int *incY)
{
    cblas_zdotc_sub(*N, X, *incX, Y, *incY, dotu);
}
void WRAP_F77(zdotu)(void *dotu, const int *N, const void *X, const int *incX,
		     const void *Y, const int *incY)
{
    cblas_zdotu_sub(*N, X, *incX, Y, *incY, dotu);
}

float WRAP_F77(sdot)(const int *N, const float *X, const int *incX,
		     const float *Y, const int *incY)
{
    return cblas_sdot(*N, X, *incX, Y, *incY);
}

float WRAP_F77(sasum)(const int *N, const float *X, const int *incX)
{
    return cblas_sasum(*N, X, *incX);
}

float WRAP_F77(scasum)(const int *N, const void *X, const int *incX)
{
    return cblas_scasum(*N, X, *incX);
}

float WRAP_F77(snrm2)(const int *N, const float *X, const int *incX)
{
    return cblas_snrm2(*N, X, *incX);
}

float WRAP_F77(scnrm2)(const int *N, const void *X, const int *incX)
{
    return cblas_scnrm2(*N, X, *incX);
}

