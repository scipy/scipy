#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H
#include "scipy_blas_defines.h"

#ifdef __cplusplus
extern "C" {
#endif

double BLAS_FUNC(dlamch)(char *cmach);

void BLAS_FUNC(dstevr)(char *jobz, char *range, CBLAS_INT *n, double *d,
                       double *e, double *vl, double *vu, CBLAS_INT *il,
                       CBLAS_INT *iu, double *abstol, CBLAS_INT *m, double *w,
                       double *z, CBLAS_INT *ldz, CBLAS_INT *issupz,
                       double *work, CBLAS_INT *lwork, CBLAS_INT *iwork,
                       CBLAS_INT *liwork, CBLAS_INT *info);

#ifdef __cplusplus
}
#endif

#endif
