/*
 * LAPACK declarations.
 */
#ifndef _SCIPY_COMMON_ARRAY_UTILS_H
#define _SCIPY_COMMON_ARRAY_UTILS_H
#include "Python.h"
#include <tuple>
#include "numpy/npy_math.h"
#include "npy_cblas.h"

/*
 * declare LAPACK prototypes
 */

extern "C" {

/* ?GETRF */
CBLAS_INT
BLAS_FUNC(sgetrf)(CBLAS_INT *m, CBLAS_INT *n, float a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(dgetrf)(CBLAS_INT *m, CBLAS_INT *n, double a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(cgetrf)(CBLAS_INT *m, CBLAS_INT *n, npy_complex64 a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(zgetrf)(CBLAS_INT *m, CBLAS_INT *n, npy_complex128 a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);


/* ?GETRI */
CBLAS_INT
BLAS_FUNC(sgetri)(CBLAS_INT *n, float a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  float work[], CBLAS_INT *lwork, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(dgetri)(CBLAS_INT *n, double a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  double work[], CBLAS_INT *lwork, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(cgetri)(CBLAS_INT *n, npy_complex64 a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  npy_complex64 work[], CBLAS_INT *lwork, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(zgetri)(CBLAS_INT *n, npy_complex128 a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  npy_complex128 work[], CBLAS_INT *lwork, CBLAS_INT *info
);


/* ?GETRS */
void BLAS_FUNC(sgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, CBLAS_INT *ipiv, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(cgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(zgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);


/* ?GECON */
void BLAS_FUNC(sgecon)(char* norm, CBLAS_INT* n, float* a,       CBLAS_INT* lda, float* anorm,  float* rcond,  float* work,       CBLAS_INT* iwork, CBLAS_INT* info);
void BLAS_FUNC(dgecon)(char* norm, CBLAS_INT* n, double* a,      CBLAS_INT* lda, double* anorm, double* rcond, double* work,      CBLAS_INT* iwork, CBLAS_INT* info);
void BLAS_FUNC(cgecon)(char* norm, CBLAS_INT* n, npy_complex64* a,  CBLAS_INT* lda, float* anorm,  float* rcond,  npy_complex64* work,  float* rwork,     CBLAS_INT* info);
void BLAS_FUNC(zgecon)(char* norm, CBLAS_INT* n, npy_complex128* a, CBLAS_INT* lda, double* anorm, double* rcond, npy_complex128* work, double* rwork,    CBLAS_INT* info);


/* ?TRTRI */
void BLAS_FUNC(strtri)(char *uplo, char *diag, CBLAS_INT *n, float* a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(dtrtri)(char *uplo, char *diag, CBLAS_INT *n, double* a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(ctrtri)(char *uplo, char *diag, CBLAS_INT *n, npy_complex64* a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(ztrtri)(char *uplo, char *diag, CBLAS_INT *n, npy_complex128* a, CBLAS_INT *lda, CBLAS_INT *info);

/* ?TRCON */
void BLAS_FUNC(strcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, float *a, CBLAS_INT *lda, float *rcond, float *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(dtrcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *rcond, double *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(ctrcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, float *rcond, npy_complex64 *work, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(ztrcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, double *rcond, npy_complex128 *work, double *rwork, CBLAS_INT *info);

/* ?TRTRS */
void BLAS_FUNC(strtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dtrtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(ctrtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(ztrtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);

/* ?POTRF */
void BLAS_FUNC(spotrf)(char *uplo, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(dpotrf)(char *uplo, CBLAS_INT *n, double *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(cpotrf)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(zpotrf)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *info);

/* ?POTRI */
void BLAS_FUNC(spotri)(char *uplo, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(dpotri)(char *uplo, CBLAS_INT *n, double *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(cpotri)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(zpotri)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *info);

/* ?POCON */
void BLAS_FUNC(spocon)(char *uplo, CBLAS_INT *n, float* a, CBLAS_INT *lda, float *anorm, float *rcond, float* work, CBLAS_INT* iwork, CBLAS_INT *info);
void BLAS_FUNC(dpocon)(char *uplo, CBLAS_INT *n, double* a, CBLAS_INT *lda, double *anorm, double *rcond, double* work, CBLAS_INT* iwork, CBLAS_INT *info);
void BLAS_FUNC(cpocon)(char *uplo, CBLAS_INT *n, npy_complex64* a, CBLAS_INT *lda, float *anorm, float *rcond, npy_complex64* work, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(zpocon)(char *uplo, CBLAS_INT *n, npy_complex128* a, CBLAS_INT *lda, double *anorm, double *rcond, npy_complex128* work, double *rwork, CBLAS_INT *info);

/* ?POTRS*/
void BLAS_FUNC(spotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dpotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(cpotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(zpotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);

/* ?SYTRF*/
void BLAS_FUNC(ssytrf)(char *uplo, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dsytrf)(char *uplo, CBLAS_INT *n, double *a, CBLAS_INT *lda, CBLAS_INT *ipiv, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(csytrf)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(zsytrf)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *work, CBLAS_INT *lwork, CBLAS_INT *info);

/* ?SYTRI */
void BLAS_FUNC(ssytri)(char *uplo, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *work, CBLAS_INT *info);
void BLAS_FUNC(dsytri)(char *uplo, CBLAS_INT *n, double *a, CBLAS_INT *lda, CBLAS_INT *ipiv, double *work, CBLAS_INT *info);
void BLAS_FUNC(csytri)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *work, CBLAS_INT *info);
void BLAS_FUNC(zsytri)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *work, CBLAS_INT *info);

/* ?SYCON*/
void BLAS_FUNC(ssycon)(char *uplo, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *anorm, float *rcond,  float *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(dsycon)(char *uplo, CBLAS_INT *n, double *a,CBLAS_INT *lda, CBLAS_INT *ipiv, double *anorm, double *rcond, double *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(csycon)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *anorm, float *rcond, npy_complex64 *work, CBLAS_INT *info);
void BLAS_FUNC(zsycon)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, double *anorm, double *rcond, npy_complex128 *work, CBLAS_INT *info);

/* ?HETRF */
void BLAS_FUNC(chetrf)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(zhetrf)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *work, CBLAS_INT *lwork, CBLAS_INT *info);

/* ?HETRI */
void BLAS_FUNC(chetri)(char *uplo, CBLAS_INT *n, npy_complex64 *a,  CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *work, CBLAS_INT *info);
void BLAS_FUNC(zhetri)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *work, CBLAS_INT *info);

/* ?SYTRS*/
void BLAS_FUNC(ssytrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dsytrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, CBLAS_INT *ipiv, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(csytrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(zsytrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);

/* ?HECON*/
void BLAS_FUNC(checon)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *anorm, float *rcond, npy_complex64 *work, CBLAS_INT *info);
void BLAS_FUNC(zhecon)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, double *anorm, double *rcond, npy_complex128 *work, CBLAS_INT *info);

/* ?HETRS*/
void BLAS_FUNC(chetrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(zhetrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);


/* ?GTTRF */
void BLAS_FUNC(sgttrf)(CBLAS_INT *n, float *dl, float *d, float *du, float *du2, CBLAS_INT *ipiv, CBLAS_INT *info);
void BLAS_FUNC(dgttrf)(CBLAS_INT *n, double *dl, double *d, double *du, double *du2, CBLAS_INT *ipiv, CBLAS_INT *info);
void BLAS_FUNC(cgttrf)(CBLAS_INT *n, npy_complex64 *dl, npy_complex64 *d, npy_complex64 *du, npy_complex64 *du2, CBLAS_INT *ipiv, CBLAS_INT *info);
void BLAS_FUNC(zgttrf)(CBLAS_INT *n, npy_complex128 *dl, npy_complex128 *d, npy_complex128 *du, npy_complex128 *du2, CBLAS_INT *ipiv, CBLAS_INT *info);

/* ?GTTRS */
void BLAS_FUNC(sgttrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, float *dl, float *d, float *du, float *du2, CBLAS_INT *ipiv, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dgttrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, double *dl, double *d, double *du, double *du2, CBLAS_INT *ipiv, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(cgttrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *dl, npy_complex64 *d, npy_complex64 *du, npy_complex64 *du2, CBLAS_INT *ipiv, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(zgttrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *dl, npy_complex128 *d, npy_complex128 *du, npy_complex128 *du2, CBLAS_INT *ipiv, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);


/* ?GTCON */
void BLAS_FUNC(sgtcon)(char *norm, CBLAS_INT *n, float *dl, float *d, float *du, float *du2, CBLAS_INT *ipiv, float *anorm, float *rcond, float *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(dgtcon)(char *norm, CBLAS_INT *n, double *dl, double *d, double *du, double *du2, CBLAS_INT *ipiv, double *anorm, double *rcond, double *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(cgtcon)(char *norm, CBLAS_INT *n, npy_complex64 *dl, npy_complex64 *d, npy_complex64 *du, npy_complex64 *du2, CBLAS_INT *ipiv, float *anorm, float *rcond, npy_complex64 *work, CBLAS_INT *info);
void BLAS_FUNC(zgtcon)(char *norm, CBLAS_INT *n, npy_complex128 *dl, npy_complex128 *d, npy_complex128 *du, npy_complex128 *du2, CBLAS_INT *ipiv, double *anorm, double *rcond, npy_complex128 *work, CBLAS_INT *info);


/* ?GESVD*/
void BLAS_FUNC(sgesvd)(char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, float *a, CBLAS_INT *lda, float *s, float *u, CBLAS_INT *ldu, float *vt, CBLAS_INT *ldvt, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dgesvd)(char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *s, double *u, CBLAS_INT *ldu, double *vt, CBLAS_INT *ldvt, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(cgesvd)(char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, float *s, npy_complex64 *u, CBLAS_INT *ldu, npy_complex64 *vt, CBLAS_INT *ldvt, npy_complex64 *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(zgesvd)(char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, double *s, npy_complex128 *u, CBLAS_INT *ldu, npy_complex128 *vt, CBLAS_INT *ldvt, npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *info);


/* ?GESDD*/
void BLAS_FUNC(sgesdd)(char *jobz, CBLAS_INT *m, CBLAS_INT *n, float *a, CBLAS_INT *lda, float *s, float *u, CBLAS_INT *ldu, float *vt, CBLAS_INT *ldvt, float *work, CBLAS_INT *lwork, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(dgesdd)(char *jobz, CBLAS_INT *m, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *s, double *u, CBLAS_INT *ldu, double *vt, CBLAS_INT *ldvt, double *work, CBLAS_INT *lwork, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(cgesdd)(char *jobz, CBLAS_INT *m, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, float *s, npy_complex64 *u, CBLAS_INT *ldu, npy_complex64 *vt, CBLAS_INT *ldvt, npy_complex64 *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(zgesdd)(char *jobz, CBLAS_INT *m, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, double *s, npy_complex128 *u, CBLAS_INT *ldu, npy_complex128 *vt, CBLAS_INT *ldvt, npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *iwork, CBLAS_INT *info);


/* ?GEQRF */
void BLAS_FUNC(sgeqrf)(CBLAS_INT *m, CBLAS_INT *n, float *a, CBLAS_INT *lda, float *tau, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dgeqrf)(CBLAS_INT *m, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *tau, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(cgeqrf)(CBLAS_INT *m, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *tau, npy_complex64 *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(zgeqrf)(CBLAS_INT *m, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *tau, npy_complex128 *work, CBLAS_INT *lwork, CBLAS_INT *info);


/* ?GEQP3 */
void BLAS_FUNC(sgeqp3)(CBLAS_INT *m, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *jpvt, float *tau, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dgeqp3)(CBLAS_INT *m, CBLAS_INT *n, double *a, CBLAS_INT *lda, CBLAS_INT *jpvt, double *tau, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(cgeqp3)(CBLAS_INT *m, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *jpvt, npy_complex64 *tau, npy_complex64 *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(zgeqp3)(CBLAS_INT *m, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *jpvt, npy_complex128 *tau, npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *info);


/* ?ORGQR, ?UNGQR */
void BLAS_FUNC(sorgqr)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *k, float *a, CBLAS_INT *lda, float *tau, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dorgqr)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *k, double *a, CBLAS_INT *lda, double *tau, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(cungqr)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *k, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *tau, npy_complex64 *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(zungqr)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *k, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *tau, npy_complex128 *work, CBLAS_INT *lwork, CBLAS_INT *info);



/* ?GELSS*/
void BLAS_FUNC(sgelss)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, float *s, float *rcond, CBLAS_INT *rank, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dgelss)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, double *s, double *rcond, CBLAS_INT *rank, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(cgelss)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a,  CBLAS_INT *lda, npy_complex64 *b,  CBLAS_INT *ldb, float *s,  float *rcond,  CBLAS_INT *rank, npy_complex64 *work,  CBLAS_INT *lwork, float *rwork,  CBLAS_INT *info);
void BLAS_FUNC(zgelss)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *b, CBLAS_INT *ldb, double *s, double *rcond, CBLAS_INT *rank, npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *info);


/* ?GELSD*/
void BLAS_FUNC(sgelsd)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, float *s, float *rcond, CBLAS_INT *rank, float *work, CBLAS_INT *lwork, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(dgelsd)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, double *s, double *rcond, CBLAS_INT *rank, double *work, CBLAS_INT *lwork, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(cgelsd)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *b, CBLAS_INT *ldb, float *s, float *rcond, CBLAS_INT *rank, npy_complex64 *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(zgelsd)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *b, CBLAS_INT *ldb, double *s, double *rcond, CBLAS_INT *rank, npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *iwork, CBLAS_INT *info);


/* ?GELSY*/
void BLAS_FUNC(sgelsy)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, float *rcond, CBLAS_INT *rank, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dgelsy)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, double *rcond, CBLAS_INT *rank, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(cgelsy)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, float *rcond, CBLAS_INT *rank, npy_complex64 *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(zgelsy)(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, double *rcond, CBLAS_INT *rank, npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *info);


/* ?GEEV, non-symmetric eigenvalues */

typedef npy_complex64 c64_t;
typedef npy_complex128 c128_t;

void BLAS_FUNC(sgeev)(char *jobvl, char *jobvr, CBLAS_INT *n, float *a,  CBLAS_INT *lda, float *wr,  float *wi,  float *vl,  CBLAS_INT *ldvl, float *vr,  CBLAS_INT *ldvr, float *work,  CBLAS_INT *lwork,                CBLAS_INT *info);
void BLAS_FUNC(dgeev)(char *jobvl, char *jobvr, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *wr, double *wi, double *vl, CBLAS_INT *ldvl, double *vr, CBLAS_INT *ldvr, double *work, CBLAS_INT *lwork,                CBLAS_INT *info);
void BLAS_FUNC(cgeev)(char *jobvl, char *jobvr, CBLAS_INT *n, c64_t *a,  CBLAS_INT *lda, c64_t *w,               c64_t *vl,  CBLAS_INT *ldvl, c64_t *vr,  CBLAS_INT *ldvr, c64_t *work,  CBLAS_INT *lwork, float *rwork,  CBLAS_INT *info);
void BLAS_FUNC(zgeev)(char *jobvl, char *jobvr, CBLAS_INT *n, c128_t *a, CBLAS_INT *lda, c128_t *w,              c128_t *vl, CBLAS_INT *ldvl, c128_t *vr, CBLAS_INT *ldvr, c128_t *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *info);


/* ?GGEV, generalized eigenvalue problem */
void BLAS_FUNC(sggev)(char *jobvl, char *jobvr, CBLAS_INT *n, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, float *alphar, float *alphai, float *beta, float *vl, CBLAS_INT *ldvl, float *vr, CBLAS_INT *ldvr, float *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(dggev)(char *jobvl, char *jobvr, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, double *alphar, double *alphai, double *beta, double *vl, CBLAS_INT *ldvl, double *vr, CBLAS_INT *ldvr, double *work, CBLAS_INT *lwork, CBLAS_INT *info);
void BLAS_FUNC(cggev)(char *jobvl, char *jobvr, CBLAS_INT *n, c64_t *a, CBLAS_INT *lda, c64_t *b, CBLAS_INT *ldb, c64_t *alpha, c64_t *beta, c64_t *vl, CBLAS_INT *ldvl, c64_t *vr, CBLAS_INT *ldvr, c64_t *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(zggev)(char *jobvl, char *jobvr, CBLAS_INT *n, c128_t *a, CBLAS_INT *lda, c128_t *b, CBLAS_INT *ldb, c128_t *alpha, c128_t *beta, c128_t *vl, CBLAS_INT *ldvl, c128_t *vr, CBLAS_INT *ldvr, c128_t *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *info);


} // extern "C"


/*
 * Generate type overloads, to map from C array types (float, double, npy_complex64, npy_complex128)
 * to LAPACK prefixes, "sdcz".
 */
#define GEN_GETRF(PREFIX, TYPE) \
inline void \
call_getrf(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getrf)(m, n, a, lda, ipiv, info); \
};

GEN_GETRF(s,float)
GEN_GETRF(d,double)
GEN_GETRF(c,npy_complex64)
GEN_GETRF(z,npy_complex128)


#define GEN_GETRS(PREFIX, TYPE) \
inline void \
call_getrs(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getrs)(trans, n, nrhs, a, lda, ipiv, b, ldb, info); \
};

GEN_GETRS(s,float)
GEN_GETRS(d,double)
GEN_GETRS(c,npy_complex64)
GEN_GETRS(z,npy_complex128)


#define GEN_GETRI(PREFIX, TYPE) \
inline void \
call_getri(CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getri)(n, a, lda, ipiv, work, lwork, info); \
};

GEN_GETRI(s,float)
GEN_GETRI(d,double)
GEN_GETRI(c,npy_complex64)
GEN_GETRI(z,npy_complex128)


// NB: iwork for real arrays or rwork for complex arrays
#define GEN_GECON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_gecon(char* norm, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## gecon)(norm, n, a, lda, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_GECON(s, float, float, CBLAS_INT)
GEN_GECON(d, double, double, CBLAS_INT)
GEN_GECON(c, npy_complex64, float, float)
GEN_GECON(z, npy_complex128, double, double)


#define GEN_TRTRI(PREFIX, TYPE) \
inline void \
call_trtri(char* uplo, char *diag, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## trtri)(uplo, diag, n, a, lda, info); \
};

GEN_TRTRI(s, float)
GEN_TRTRI(d, double)
GEN_TRTRI(c, npy_complex64)
GEN_TRTRI(z, npy_complex128)


#define GEN_TRCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_trcon(char* norm, char *uplo, char *diag, CBLAS_INT *n, CTYPE *a, CBLAS_INT *lda, RTYPE *rcond, CTYPE *work, void *irwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## trcon)(norm, uplo, diag, n, a, lda, rcond, work, (WTYPE *)irwork, info); \
};

GEN_TRCON(s, float, float, CBLAS_INT)
GEN_TRCON(d, double, double, CBLAS_INT)
GEN_TRCON(c, npy_complex64, float, float)
GEN_TRCON(z, npy_complex128, double, double)


#define GEN_TRTRS(PREFIX, TYPE) \
inline void \
call_trtrs(char* uplo, char *trans, char *diag, CBLAS_INT* n, CBLAS_INT* nrhs, TYPE* a, CBLAS_INT* lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## trtrs)(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info); \
};

GEN_TRTRS(s, float)
GEN_TRTRS(d, double)
GEN_TRTRS(c, npy_complex64)
GEN_TRTRS(z, npy_complex128)


#define GEN_POTRF(PREFIX, TYPE) \
inline void \
call_potrf(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potrf)(uplo, n, a, lda, info); \
};

GEN_POTRF(s, float)
GEN_POTRF(d, double)
GEN_POTRF(c, npy_complex64)
GEN_POTRF(z, npy_complex128)


#define GEN_POTRI(PREFIX, TYPE) \
inline void \
call_potri(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potri)(uplo, n, a, lda, info); \
};

GEN_POTRI(s, float)
GEN_POTRI(d, double)
GEN_POTRI(c, npy_complex64)
GEN_POTRI(z, npy_complex128)


// NB: iwork for real arrays or rwork for complex arrays
#define GEN_POCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_pocon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## pocon)(uplo, n, a, lda, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_POCON(s, float, float, CBLAS_INT)
GEN_POCON(d, double, double, CBLAS_INT)
GEN_POCON(c, npy_complex64, float, float)
GEN_POCON(z, npy_complex128, double, double)


#define GEN_POTRS(PREFIX, TYPE) \
inline void \
call_potrs(char* uplo, CBLAS_INT* n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potrs)(uplo, n, nrhs, a, lda, b, ldb, info); \
};

GEN_POTRS(s, float)
GEN_POTRS(d, double)
GEN_POTRS(c, npy_complex64)
GEN_POTRS(z, npy_complex128)


#define GEN_SYTRF(PREFIX, TYPE) \
inline void \
call_sytrf(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sytrf)(uplo, n, a, lda, ipiv, work, lwork, info); \
};

GEN_SYTRF(s, float)
GEN_SYTRF(d, double)
GEN_SYTRF(c, npy_complex64)
GEN_SYTRF(z, npy_complex128)


// dispatch to sSYtrf for "float hermitian"
#define GEN_HETRF(PREFIX, L_PREFIX, TYPE) \
inline void \
call_hetrf(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## L_PREFIX ## trf)(uplo, n, a, lda, ipiv, work, lwork, info); \
};

GEN_HETRF(s, sy, float)
GEN_HETRF(d, sy, double)
GEN_HETRF(c, he, npy_complex64)
GEN_HETRF(z, he, npy_complex128)


#define GEN_SYTRI(PREFIX, TYPE) \
inline void \
call_sytri(char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## sytri)(uplo, n, a, lda, ipiv, work, info); \
};

GEN_SYTRI(s, float)
GEN_SYTRI(d, double)
GEN_SYTRI(c, npy_complex64)
GEN_SYTRI(z, npy_complex128)


// dispatch to sSYtri for "float hermitian"
#define GEN_HETRI(PREFIX, L_PREFIX, TYPE) \
inline void \
call_hetri(char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## L_PREFIX ## tri)(uplo, n, a, lda, ipiv, work, info); \
};

GEN_HETRI(s, sy, float)
GEN_HETRI(d, sy, double)
GEN_HETRI(c, he, npy_complex64)
GEN_HETRI(z, he, npy_complex128)


// NB: iwork for real arrays only, no rwork for complex routines (10 arguments for s- d- variants; 9 arguments for c- and z- variants)
#define GEN_SYCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_sycon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sycon)(uplo, n, a, lda, ipiv, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_SYCON(s, float, float, CBLAS_INT)
GEN_SYCON(d, double, double, CBLAS_INT)

#define GEN_SYCON_CZ(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_sycon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sycon)(uplo, n, a, lda, ipiv, anorm, rcond, work, info); \
};

GEN_SYCON_CZ(c, npy_complex64, float, float)
GEN_SYCON_CZ(z, npy_complex128, double, double)


// dispatch to sSYcon for "float hermitian"
#define GEN_HECON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_hecon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sycon)(uplo, n, a, lda, ipiv, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_HECON(s, float, float, CBLAS_INT)
GEN_HECON(d, double, double, CBLAS_INT)

#define GEN_HECON_CZ(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_hecon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## hecon)(uplo, n, a, lda, ipiv, anorm, rcond, work, info); \
};

GEN_HECON_CZ(c, npy_complex64, float, float)
GEN_HECON_CZ(z, npy_complex128, double, double)


#define GEN_SYTRS(PREFIX, TYPE) \
inline void \
call_sytrs(char* uplo, CBLAS_INT* n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sytrs)(uplo, n, nrhs, a, lda, ipiv, b, ldb, info); \
};

GEN_SYTRS(s, float)
GEN_SYTRS(d, double)
GEN_SYTRS(c, npy_complex64)
GEN_SYTRS(z, npy_complex128)


// dispatch to sSYtrs for "float hermitian"
#define GEN_HETRS(PREFIX, L_PREFIX, TYPE) \
inline void \
call_hetrs(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## L_PREFIX ## trs)(uplo, n, nrhs, a, lda, ipiv, b, ldb, info); \
};

GEN_HETRS(s, sy, float)
GEN_HETRS(d, sy, double)
GEN_HETRS(c, he, npy_complex64)
GEN_HETRS(z, he, npy_complex128)


#define GEN_GTTRF(PREFIX, TYPE) \
inline void \
call_gttrf(CBLAS_INT *n, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gttrf)(n, dl, d, du, du2, ipiv, info); \
};

GEN_GTTRF(s, float)
GEN_GTTRF(d, double)
GEN_GTTRF(c, npy_complex64)
GEN_GTTRF(z, npy_complex128)


#define GEN_GTTRS(PREFIX, TYPE) \
inline void \
call_gttrs(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gttrs)(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info); \
};

GEN_GTTRS(s, float)
GEN_GTTRS(d, double)
GEN_GTTRS(c, npy_complex64)
GEN_GTTRS(z, npy_complex128)


#define GEN_GTCON(PREFIX, TYPE) \
inline void \
call_gtcon(char *norm, CBLAS_INT *n, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, TYPE *anorm, TYPE *rcond, TYPE *work, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gtcon)(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info); \
};

GEN_GTCON(s, float)
GEN_GTCON(d, double)


// NB: `iwork` is not used for c- and z- variants of ?gtcon
#define GEN_GTCON_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gtcon(char *norm, CBLAS_INT *n, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, RTYPE *anorm, RTYPE *rcond, TYPE *work, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gtcon)(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info); \
};

GEN_GTCON_CZ(c, npy_complex64, float)
GEN_GTCON_CZ(z, npy_complex128, double)



/*
 * ?GESVD wrappers.
 *
 * We need to wrap over:
 *   - four type variants, s-, d-, c-, and zgesvd;
 *   - complex variants, c- and z-, receive the `rwork` argument, while s- and d- variants do not.
 * Thus,
 *   - `call_gesvd` has four overloads;
 *   - all variants receive the `rwork` argument; c- and z- variants forward it to LAPACK,
 *     and s- and d- variants swallow it.
 */
inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, float *a, CBLAS_INT *lda,
    float *s, float *u, CBLAS_INT *ldu, float *vt, CBLAS_INT *ldvt, float *work, CBLAS_INT *lwork,
    float *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(sgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
};

inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, double *a, CBLAS_INT *lda,
    double *s, double *u, CBLAS_INT *ldu, double *vt, CBLAS_INT *ldvt, double *work, CBLAS_INT *lwork,
    double *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(dgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
};

inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda,
    float *s, npy_complex64 *u, CBLAS_INT *ldu, npy_complex64 *vt, CBLAS_INT *ldvt,
    npy_complex64 *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(cgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
};

inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda,
    double *s, npy_complex128 *u, CBLAS_INT *ldu, npy_complex128 *vt, CBLAS_INT *ldvt,
    npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(zgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
};



/*
 * ?GESDD wrappers.
 *
 * The logic is similar to ?gesdd:
 *   - we overload for four type variants, s-, d-, c-, and z-;
 *   - we forward `rwork` to c- and z- LAPACK functions, and swallow it for s- and d-;
 *
 */

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, float *a, CBLAS_INT *lda, float *s, float *u, CBLAS_INT *ldu,
    float *vt, CBLAS_INT *ldvt, float *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(sgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
};

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *s, double *u, CBLAS_INT *ldu,
    double *vt, CBLAS_INT *ldvt, double *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(dgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
};

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, float *s, npy_complex64 *u, CBLAS_INT *ldu,
    npy_complex64 *vt, CBLAS_INT *ldvt, npy_complex64 *work, CBLAS_INT *lwork, float *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(cgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
};

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, double *s, npy_complex128 *u, CBLAS_INT *ldu,
    npy_complex128 *vt, CBLAS_INT *ldvt, npy_complex128 *work, CBLAS_INT *lwork, double *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(zgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
};




#define GEN_GEQRF(PREFIX, TYPE) \
inline void \
geqrf(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *tau, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## geqrf)(m, n, a, lda, tau, work, lwork, info); \
};

GEN_GEQRF(s, float)
GEN_GEQRF(d, double)
GEN_GEQRF(c, npy_complex64)
GEN_GEQRF(z, npy_complex128)



// NB: wrap {s-,d-}orgqr for reals and {c-,z-}ungqr for complex
#define GEN_OR_UN_GQR(PREFIX, TYPE) \
inline void \
or_un_gqr(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *k, TYPE *a, CBLAS_INT *lda, TYPE *tau, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gqr)(m, n, k, a, lda, tau, work, lwork, info); \
};

GEN_OR_UN_GQR(sor, float)
GEN_OR_UN_GQR(dor, double)
GEN_OR_UN_GQR(cun, npy_complex64)
GEN_OR_UN_GQR(zun, npy_complex128)



// NB: s- and d- variants ignore the rwork argument (because LAPACK routines do not have it
#define GEN_GELSS_SD(PREFIX, TYPE, RTYPE) \
inline void \
call_gelss(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelss)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info); \
};

GEN_GELSS_SD(s, float, float)
GEN_GELSS_SD(d, double, double)


#define GEN_GELSS_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gelss(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelss)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info); \
};

GEN_GELSS_CZ(c, npy_complex64, float)
GEN_GELSS_CZ(z, npy_complex128, double)


// NB: s- and d- variants ignore the rwork argument (because LAPACK routines do not have it
#define GEN_GELSD_SD(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsd(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsd)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info); \
};

GEN_GELSD_SD(s, float, float)
GEN_GELSD_SD(d, double, double)


#define GEN_GELSD_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsd(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsd)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info); \
};

GEN_GELSD_CZ(c, npy_complex64, float)
GEN_GELSD_CZ(z, npy_complex128, double)


// NB: s- and d- variants ignore the rwork argument (because LAPACK routines do not have it
#define GEN_GELSY_SD(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsy(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsy)(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info); \
};

GEN_GELSY_SD(s, float, float)
GEN_GELSY_SD(d, double, double)


#define GEN_GELSY_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsy(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsy)(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info); \
};

GEN_GELSY_CZ(c, npy_complex64, float)
GEN_GELSY_CZ(z, npy_complex128, double)


/*
 * ?GEEV wrappers.
 *
 * We need to wrap over:
 *   - four type variants, s-, d-, c-, and zgeev;
 *   - complex variants, c- and z-, receive the `rwork` argument, while s- and d- variants do not.
 *   - s- and d- variants return real and imaginary parts of eigenvalues separately, in *wr and *wi arrays
 *     c- and z- variants return a single complex array, *w, instead
 * Thus,
 *   - `call_geev` has four overloads;
 *   - all variants receive the `rwork` argument; c- and z- variants forward it to LAPACK,
 *     and s- and d- variants swallow it.
 *   - all variants have *wr and *wi arguments, both of the same type as *a
 *     (real for real *a, complex for complex *a);
 *     real-valued overloads, s- and d-, only fill *wr and ignore the *wi argument.
 */
#define GEN_GEEV_SD(PREFIX, TYPE) \
inline void \
call_geev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *wr, TYPE *wi, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork,  TYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## geev)(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info); \
};

GEN_GEEV_SD(s, float)
GEN_GEEV_SD(d, double)

#define GEN_GEEV_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_geev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *wr, TYPE *wi, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    /* ignore wi */ \
    BLAS_FUNC(PREFIX ## geev)(jobvl, jobvr, n, a, lda, wr, vl, ldvl, vr, ldvr, work, lwork, rwork, info); \
};

GEN_GEEV_CZ(c, npy_complex64, float)
GEN_GEEV_CZ(z, npy_complex128, double)


/*
 * Wrappers for ?GGEV
 *
 * The design is similar to that of ?GEEV wrappers: all overloads receive *rwork and *alphar, *alphai,
 */
#define GEN_GGEV_SD(PREFIX, TYPE) \
inline void \
call_ggev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, TYPE *alphar, TYPE *alphai, TYPE *beta, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork, TYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## ggev)(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info); \
};

GEN_GGEV_SD(s, float)
GEN_GGEV_SD(d, double)


#define GEN_GGEV_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_ggev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, TYPE *alphar, TYPE *alphai, TYPE *beta, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## ggev)(jobvl, jobvr, n, a, lda, b, ldb, alphar, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info); \
};

GEN_GGEV_CZ(c, npy_complex64, float)
GEN_GGEV_CZ(z, npy_complex128, double)


// Structure tags; python side maps assume_a strings to these values
enum St : Py_ssize_t
{
    NONE = -1,
    GENERAL = 0,
    DIAGONAL = 11,
    TRIDIAGONAL = 31,
    UPPER_TRIANGULAR = 21,
    LOWER_TRIANGULAR = 22,
    POS_DEF = 101,
    SYM = 201,
    HER = 211
};


/*
 * Rich return object
 */
struct SliceStatus {
    Py_ssize_t slice_num;
    Py_ssize_t structure;
    int is_singular;
    int is_ill_conditioned;
    double rcond;
    Py_ssize_t lapack_info;
};


void init_status(SliceStatus& slice_status, npy_intp idx, St slice_structure) {
    slice_status.slice_num = idx;
    slice_status.structure = (Py_ssize_t)slice_structure;
    slice_status.is_singular = 0;
    slice_status.is_ill_conditioned = 0;
    slice_status.rcond = 0;
    slice_status.lapack_info = 0;
}


typedef std::vector<SliceStatus> SliceStatusVec;


/*
 * When looping over slices, each `solve_slice_XYZ` return the `status` struct,
 * which records one of three possible outcomes:
 *
 *  - the slice is singular or LAPACK returned non-zero `info`
 *  - the slice is detected to be ill-conditioned
 *  - all is well
 *
 * In the first two cases, we record the `status` in `vec_status`.
 * For non-recoverable errors (singularity of non-zero info), we terminate the loop over the slices.
 */
int
_detect_problems(const SliceStatus& slice_status, SliceStatusVec& vec_status) {
    if ((slice_status.lapack_info < 0) || (slice_status.is_singular)) {
        vec_status.push_back(slice_status);
        return 1;
    }
    else if (slice_status.is_ill_conditioned) {
        vec_status.push_back(slice_status);
    }
    return 0;
}


/*
 * lwork defensive handler :
 *  cf https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/lapack.py#L1004
 *
 *  Round floating-point lwork returned by lapack to integer.
 *
 *  Several LAPACK routines compute optimal values for LWORK, which
 *  they return in a floating-point variable. However, for large
 *  values of LWORK, single-precision floating point is not sufficient
 *  to hold the exact value --- some LAPACK versions (<= 3.5.0 at
 *  least) truncate the returned integer to single precision and in
 *  some cases this can be smaller than the required value.
 *
 *  The fudge_factor comes from
 *  https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/_basic.py#L1154
 *
 *  A fudge factor of 1% was added in commit https://github.com/scipy/scipy/commit/dfb543c147c
 *  to avoid a "curious segfault with 500x500 matrices and OpenBLAS".
 */
template<typename T>
CBLAS_INT _calc_lwork(T _lwrk, double fudge_factor=1.0) {
    using real_type = typename type_traits<T>::real_type;

    real_type value = real_part(_lwrk) * fudge_factor;
    if((std::is_same<real_type, float>::value) ||
       (std::is_same<real_type, npy_complex64>::value)
    ) {
        // Single-precision routine -- take next fp value to work
        // around possible truncation in LAPACK code
        value = std::nextafter(value, std::numeric_limits<real_type>::infinity());
    }

    CBLAS_INT lwork;
    if ((value < 0) || !(value <= (real_type)std::numeric_limits<CBLAS_INT>::max())) {
        // Too large lwork required - Computation cannot be performed with standard LAPACK
        lwork = -1;
    }
    else {
        lwork = value > 0 ? (CBLAS_INT)value : 1;
    }
    return lwork;
}


/*
 * Given an array of ndim >= 2, compute a pointer to the start of the 2D "slice" idx
 */
template<typename T>
T* compute_slice_ptr(npy_intp idx, T *Am_data, npy_intp ndim, npy_intp *shape, npy_intp *strides) {
    npy_intp offset = 0;
    npy_intp temp_idx = idx;
    for (int i = ndim - 3; i >= 0; i--) {
        offset += (temp_idx % shape[i]) * strides[i];
        temp_idx /= shape[i];
    }
    T* slice_ptr = (T *)(Am_data + (offset/sizeof(T)));
    return slice_ptr;
}


/*
 * Copy n-by-m slice from slice_ptr to dst.
 */
template<typename T>
void copy_slice(T* dst, const T* slice_ptr, const npy_intp n, const npy_intp m, const npy_intp s2, const npy_intp s1) {

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i * m + j] = *(slice_ptr + (i*s2/sizeof(T)) + (j*s1/sizeof(T)));
        }
    }
}


/*
 * Copy n-by-m C-order slice from slice_ptr to dst in F-order.
 *
 * `src` is n-by-m, strided
 * `dst` is ldb-by-m, F-ordered.
 *
 * The default is to have src and dst of the same size (ldb=-1 means ldb=n).
 */
template<typename T>
void copy_slice_F(T* dst, const T* slice_ptr, const npy_intp n, const npy_intp m, const npy_intp s2, const npy_intp s1, npy_intp ldb=-1) {

    if (ldb == -1) {ldb = n;}

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i + j*ldb] = *(slice_ptr + (i*s2/sizeof(T)) + (j*s1/sizeof(T)));  // == src[i*m + j]
        }
    }
}

/*
 * Copy n-by-m F-ordered `src` to C-ordered `dst`.
 *
 * `src` is ldb-by-m, F-ordered
 * `dst` is n-by-m, C-ordered
 *
 * The default is to have src and dst of the same size (ldb=-1 means ldb=n)
 */
template<typename T>
void copy_slice_F_to_C(T* dst, const T* src, const npy_intp n, const npy_intp m, npy_intp ldb=-1) {

    if (ldb == -1) {ldb = n;}

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i*m + j] = src[i + j*ldb];
        }
    }
}





/*
 * 1-norm of a matrix
 */

template<typename T>
typename type_traits<T>::real_type
norm1_(T* A, T* work, const npy_intp n)
{
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;
    value_type *pA = reinterpret_cast<value_type *>(A);

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { rwork[i] = std::abs(pA[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { rwork[j] += std::abs(pA[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
typename type_traits<T>::real_type
norm1_sym_herm_upper(T* A, T* work, const npy_intp n)
{
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;
    value_type *pA = reinterpret_cast<value_type *>(A);

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { rwork[i] = std::abs(pA[i]);
     }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) {
        // only loop over the upper triangle
        rwork[i] += std::abs(pA[i*n + i]);
        for (j = i+1; j < n; j++) {
            temp = std::abs(pA[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
typename type_traits<T>::real_type
norm1_sym_herm_lower(T* A, T* work, const npy_intp n)
{
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;
    value_type *pA = reinterpret_cast<value_type *>(A);

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    for (i = 0; i < n; i++) { rwork[i] = 0.0; }

    for (i=0; i < n; i++) {
        rwork[i] += std::abs(pA[i*n + i]);
        for (j=0; j < i; j++) {
            temp = std::abs(pA[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }

    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
typename type_traits<T>::real_type
norm1_sym_herm(char uplo, T *A, T *work, const npy_intp n) {
    // NB: transpose for the F order
    if (uplo == 'U') {return norm1_sym_herm_lower(A, work, n);}
    else if (uplo == 'L') {return norm1_sym_herm_upper(A, work, n);}
    else {throw std::runtime_error("uplo at norms");}
}


template<typename T>
typename type_traits<T>::real_type
norm1_tridiag(T* dl, T *d, T *du, T *work, const npy_intp n) {
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;

    value_type *pd = reinterpret_cast<value_type *>(d);
    value_type *pdu = reinterpret_cast<value_type *>(du);
    value_type *pdl = reinterpret_cast<value_type *>(dl);
    real_type *rwork = (real_type *)work;

    npy_intp i;
    for (i=0; i<n; i++) {
        rwork[i] = std::abs(pd[i]);
    }
    for (i=0; i<n-1; i++) {
        rwork[i] += std::abs(pdl[i]);
    }
    for (i=1; i<n-1; i++) {
        rwork[i] += std::abs(pdu[i-1]);
    }

    real_type temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


/***************************
 ***  Structure detection
 ***************************/

template<typename T>
void
bandwidth(T* data, npy_intp n, npy_intp m, npy_intp* lower_band, npy_intp* upper_band)
{
    using value_type = typename type_traits<T>::value_type;
    value_type *p_data = reinterpret_cast<value_type *>(data);
    value_type zero = value_type(0.);

    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t c = 0; c < m-1; c++)
    {
        for (Py_ssize_t r = n-1; r > c + lb; r--)
        {
            if (p_data[c*n + r] != zero) { lb = r - c; break; }
        }
        if (c + lb + 1 > m) { break; }
    }
    for (Py_ssize_t c = m-1; c > 0; c--)
    {
        for (Py_ssize_t r = 0; r < c - ub; r++)
        {
            if (p_data[c*n + r] != zero) { ub = c - r; break; }

        }
        if (c <= ub) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


template<typename T>
std::tuple<bool, bool>
is_sym_or_herm(const T *data, npy_intp n) {
    // Return a pair of (is_symmetric, is_hermitian)
    using value_type = typename type_traits<T>::value_type;
    const value_type *p_data = reinterpret_cast<const value_type *>(data);
    bool all_sym = true, all_herm = true;

    for (npy_intp i=0; i < n; i++) {
        for (npy_intp j=0; j < n; j++) {
            value_type elem1 = p_data[i*n + j];
            value_type elem2 = p_data[i + j*n];
            all_sym = all_sym && (elem1 == elem2);
            all_herm = all_herm && (elem1 == std::conj(elem2));
            if(!(all_sym || all_herm)) {
                // short-circuit : it's neither symmetric not hermitian
                return std::make_tuple(false, false);
            }
        }
    }
    return std::make_tuple(all_sym, all_herm);
}


template<typename T>
inline void
swap_cf(T* src, T* dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    T *bb = dst;
    T *aa = src;
    if ((r < 16) && (c < 16)) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf(src, dst, r2, c, n);
            swap_cf(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf(src, dst, r, c2, n);
            swap_cf(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
}


/*
 * Common matrices
 */

// fill np.triu(a) from np.tril(a) or np.tril(a) from np.triu(a)
template<typename T>
inline void
fill_other_triangle(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                data[j + i*n] = conj(data[i + j*n]);
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i+1; j++){
                data[j + i*n] = conj(data[i + j*n]);
            }
        }
    }
}

// XXX deduplicate (conj or noconj)
template<typename T>
inline void
fill_other_triangle_noconj(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                data[j + i*n] = data[i + j*n];
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i+1; j++){
                data[j + i*n] = data[i + j*n];
            }
        }
    }
}


/*
 * Helper for converting an NxN matrix to tridiagonal form:
 * extract the diagonals of `data` (a full NxN matrix) into du, d, dl
 *
 * the upper and lower subdiagonals, `dl` and `du`, are have length N-1
 * the main diagonal, `d`, has length N
 *
 */
template<typename T>
inline void
to_tridiag(const T *data, npy_intp N, T *du, T *d, T *dl) {
    for (npy_intp i=0; i<N; i++) {
        d[i] = data[i + i*N];
    }
    for (npy_intp i=0; i<N-1; i++) {
        dl[i] = data[i*N + i + 1];
    }
    for (npy_intp i=0; i<N-1; i++) {
        du[i] = data[(i+1)*N + i];
    }
}


template<typename T>
inline void
zero_other_triangle(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                data[j + i*n] = numeric_limits<T>::zero;
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i; j++){
                data[j + i*n] = numeric_limits<T>::zero;
            }
        }
    }
}


template<typename T>
inline void
nan_matrix(T * data, npy_intp n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            data[i * n + j] = numeric_limits<T>::nan;
        }
    }
}
#endif




