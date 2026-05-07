#ifndef BLASKLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H
#include "scipy_blas_defines.h"

#ifdef __cplusplus
extern "C" {
#endif

void BLAS_FUNC(dstevd)(char *jobz, CBLAS_INT *n, double *d, double *e, double *z, CBLAS_INT *ldz, double *work, CBLAS_INT *lwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info);

#ifdef __cplusplus
  }
#endif

#endif
