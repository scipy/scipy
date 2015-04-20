#include "cblas.h"

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

void WRAP_F77(acc_sdot)(const int *N, const float *X, const int *incX,
                        const float *Y, const int *incY,
                        float *result)
{
    *result = cblas_sdot(*N, X, *incX, Y, *incY);
}

void WRAP_F77(acc_sdsdot)(const int *N, const float *alpha,
                          const float *X, const int *incX,
                          const float *Y, const int *incY,
                          float *result)
{
    *result = cblas_sdsdot(*N, *alpha, X, *incX, Y, *incY);
}

void WRAP_F77(acc_sasum)(const int *N, const float *X, const int *incX,
                         float *result)
{
    *result = cblas_sasum(*N, X, *incX);
}

void WRAP_F77(acc_snrm2)(const int *N, const float *X, const int *incX,
                         float *result)
{
    *result = cblas_snrm2(*N, X, *incX);
}

void WRAP_F77(acc_scasum)(const int *N, const void *X, const int *incX,
                          float *result)
{
    *result = cblas_scasum(*N, X, *incX);
}

void WRAP_F77(acc_scnrm2)(const int *N, const void *X, const int *incX,
                          float *result)
{
    *result = cblas_scnrm2(*N, X, *incX);
}
