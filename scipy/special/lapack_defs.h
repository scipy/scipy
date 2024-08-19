/*
 * Handle different Fortran conventions.
 */

#include "npy_cblas.h"

extern void BLAS_FUNC(dstevr)(char *jobz, char *range, CBLAS_INT *n, double *d, double *e,
                              double *vl, double *vu, CBLAS_INT *il, CBLAS_INT *iu, double *abstol,
                              CBLAS_INT *m, double *w, double *z, CBLAS_INT *ldz, CBLAS_INT *isuppz,
                              double *work, CBLAS_INT *lwork, CBLAS_INT *iwork, CBLAS_INT *liwork,
                              CBLAS_INT *info, size_t jobz_len, size_t range_len);

static void c_dstevr(char *jobz, char *range, CBLAS_INT *n, double *d, double *e,
                     double *vl, double *vu, CBLAS_INT *il, CBLAS_INT *iu, double *abstol,
                     CBLAS_INT *m, double *w, double *z, CBLAS_INT *ldz, CBLAS_INT *isuppz,
                     double *work, CBLAS_INT *lwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info) {
    BLAS_FUNC(dstevr)(jobz, range, n, d, e, vl, vu, il, iu, abstol, m,
                      w, z, ldz, isuppz, work, lwork, iwork, liwork, info,
                      1, 1);
}
