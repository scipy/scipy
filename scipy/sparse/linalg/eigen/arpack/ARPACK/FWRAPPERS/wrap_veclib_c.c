#include <Accelerate/Accelerate.h>

#define WRAP_F77(a) a##_

/* This intermediate C-file is necessary also for the Fortran wrapper:
   - The cblas_* functions cannot be called from Fortran directly
     as they do not take all of their arguments as pointers
     (as is done in Fortran).
   - It would in principle be possible to make a wrapper directly between
     Fortran as for Lapack. However, apparently some versions of
     MacOSX have a bug in the FORTRAN interface to SDOT
     (). In addition, the cblas seems to be the only standard officially
     endorsed by Apple.
*/

float WRAP_F77(acc_sdot)(const int *N, const float *X, const int *incX,
			 const float *Y, const int *incY)
{
    return cblas_sdot(*N, X, *incX, Y, *incY);
}

float WRAP_F77(acc_sdsdot)(const int *N, const float *alpha,
			   const float *X, const int *incX,
			   const float *Y, const int *incY)
{
    return cblas_sdsdot(*N, *alpha, X, *incX, Y, *incY);
}

float WRAP_F77(acc_sasum)(const int *N, const float *X, const int *incX)
{
    return cblas_sasum(*N, X, *incX);
}

float WRAP_F77(acc_snrm2)(const int *N, const float *X, const int *incX)
{
    return cblas_snrm2(*N, X, *incX);
}

float WRAP_F77(acc_scasum)(const int *N, const void *X, const int *incX)
{
    return cblas_scasum(*N, X, *incX);
}

float WRAP_F77(acc_scnrm2)(const int *N, const void *X, const int *incX)
{
    return cblas_scnrm2(*N, X, *incX);
}

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

