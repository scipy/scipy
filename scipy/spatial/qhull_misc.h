/*
 * Handle different Fortran conventions and qh_new_qhull_scipy entry point.
 */
#ifndef QHULL_MISC_H_
#define QHULL_MISC_H_

#include "npy_cblas.h"

void BLAS_FUNC(dgetrf)(CBLAS_INT*, CBLAS_INT*, double*, CBLAS_INT*, CBLAS_INT*, CBLAS_INT*);
void BLAS_FUNC(dgecon)(char*, CBLAS_INT*, double*, CBLAS_INT*, double*, double*, double*, CBLAS_INT*, CBLAS_INT*, size_t);
void BLAS_FUNC(dgetrs)(char*, CBLAS_INT*, CBLAS_INT*, double*, CBLAS_INT*, CBLAS_INT*, double*, CBLAS_INT*, CBLAS_INT*, size_t);

#define qh_dgetrf(m,n,a,lda,ipiv,info) BLAS_FUNC(dgetrf)(m,n,a,lda,ipiv,info)
#define qh_dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info) BLAS_FUNC(dgecon)(norm,n,a,lda,anorm,rcond,work,iwork,info,1)
#define qh_dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info) BLAS_FUNC(dgetrs)(trans,n,nrhs,a,lda,ipiv,b,ldb,info,1)

#define qhull_misc_lib_check() QHULL_LIB_CHECK

#include "qhull_src/src/libqhull_r.h"

int qh_new_qhull_scipy(qhT *qh, int dim, int numpoints, coordT *points, boolT ismalloc,
                       char *qhull_cmd, FILE *outfile, FILE *errfile, coordT* feaspoint);

#endif /* QHULL_MISC_H_ */
