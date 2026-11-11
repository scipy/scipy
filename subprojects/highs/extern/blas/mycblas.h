#ifndef FACTORHIGHS_MY_CBLAS_H
#define FACTORHIGHS_MY_CBLAS_H

#include <stdint.h>

#if defined(HIPO_USES_APPLE_BLAS)

// Use Apple's own Accelerate header (instead of the hand-rolled Netlib-style
// declarations below) so that cblas_* calls bind to the modern, actively
// maintained Accelerate implementation instead of the legacy one (frozen at
// LAPACK 3.2.1, with long-standing known correctness quirks in its BLAS).
// ACCELERATE_NEW_LAPACK only takes effect through Apple's header itself
// (it works via $NEWLAPACK asm-label symbol renaming declared there); the
// macro alone, without including this header, has no effect - see
// https://github.com/ERGO-Code/HiGHS/issues/3160.
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>

typedef int32_t blasint;

#else

// Provide definition for cblas functions
// Based on Netlib implementation

enum CBLAS_ORDER { CblasRowMajor = 101, CblasColMajor = 102 };
enum CBLAS_TRANSPOSE {
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113
};
enum CBLAS_UPLO { CblasUpper = 121, CblasLower = 122 };
enum CBLAS_DIAG { CblasNonUnit = 131, CblasUnit = 132 };
enum CBLAS_SIDE { CblasLeft = 141, CblasRight = 142 };

#ifdef __cplusplus
extern "C" {
#endif

typedef int32_t blasint;

// level 1

void cblas_daxpy(const blasint n, const double alpha, const double* x,
                 const blasint incx, double* y, const blasint incy);
void cblas_dcopy(const blasint n, const double* x, const blasint incx,
                 double* y, const blasint incy);
void cblas_dscal(const blasint n, const double alpha, double* x,
                 const blasint incx);
void cblas_dswap(const blasint n, double* x, const blasint incx, double* y,
                 const blasint incy);

// level 2

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa, const blasint M,
                 const blasint n, const double alpha, const double* A,
                 const blasint lda, const double* x, const blasint incx,
                 const double beta, double* y, const blasint incy);

void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const blasint n, const double* ap, double* x,
                 const blasint incx);

void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const blasint n, const double* a, const blasint lda, double* x,
                 const blasint incx);

void cblas_dger(const enum CBLAS_ORDER order, const blasint m, const blasint n,
                const double alpha, const double* x, const blasint incx,
                const double* y, const blasint incy, double* A,
                const blasint lda);

// level 3

void cblas_dgemm(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_TRANSPOSE transb, const blasint m,
                 const blasint n, const blasint k, const double alpha,
                 const double* A, const blasint lda, const double* B,
                 const blasint ldb, const double beta, double* C,
                 const blasint ldc);

void cblas_dsyrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const blasint n,
                 const blasint k, const double alpha, const double* a,
                 const blasint lda, const double beta, double* C,
                 const blasint ldc);

void cblas_dtrsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
                 const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_DIAG diag, const blasint m, const blasint n,
                 const double alpha, const double* a, const blasint lda,
                 double* b, const blasint ldb);

#ifdef __cplusplus
}
#endif

#endif  // HIPO_USES_APPLE_BLAS

#ifdef __cplusplus
extern "C" {
#endif

#if defined(HIPO_USES_OPENBLAS)
void openblas_set_num_threads(int num_threads);
int openblas_get_num_threads(void);
#endif

int highs_openblas_set_num_threads(int num_threads);
#ifdef __cplusplus
}
#endif

#endif