#include "CallAndTimeBlas.h"

#include "DataCollector.h"
#include "DenseFact.h"
#include "FactorHighsSettings.h"
#include "HighsExternalApi.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

// macros to interface with CBlas
#define TRANS(x) (x) == 'N' ? CblasNoTrans : CblasTrans
#define UPLO(x) (x) == 'U' ? CblasUpper : CblasLower
#define DIAG(x) (x) == 'N' ? CblasNonUnit : CblasUnit
#define SIDE(x) (x) == 'L' ? CblasLeft : CblasRight

// level 1

void callAndTime_daxpy(Int n, double da, const double* dx, Int incx, double* dy,
                       Int incy, DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::daxpy(n, da, dx, incx, dy, incy);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_axpy);
}

void callAndTime_dcopy(Int n, const double* dx, Int incx, double* dy, Int incy,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dcopy(n, dx, incx, dy, incy);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_copy);
}

void callAndTime_dscal(Int n, const double da, double* dx, Int incx,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dscal(n, da, dx, incx);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_scal);
}

void callAndTime_dswap(Int n, double* dx, Int incx, double* dy, Int incy,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dswap(n, dx, incx, dy, incy);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_swap);
}

// level 2

void callAndTime_dgemv(char trans, Int m, Int n, double alpha, const double* A,
                       Int lda, const double* x, Int incx, double beta,
                       double* y, Int incy, DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dgemv(CblasColMajor, TRANS(trans), m, n, alpha, A, lda, x,
                           incx, beta, y, incy);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_gemv);
}

void callAndTime_dtpsv(char uplo, char trans, char diag, Int n,
                       const double* ap, double* x, Int incx,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dtpsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag),
                           n, ap, x, incx);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_tpsv);
}

void callAndTime_dtrsv(char uplo, char trans, char diag, Int n, const double* A,
                       Int lda, double* x, Int incx, DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dtrsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag),
                           n, A, lda, x, incx);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_trsv);
}

void callAndTime_dger(Int m, Int n, double alpha, const double* x, Int incx,
                      const double* y, Int incy, double* A, Int lda,
                      DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_ger);
}

// level 3

void callAndTime_dgemm(char transa, char transb, Int m, Int n, Int k,
                       double alpha, const double* A, Int lda, const double* B,
                       Int ldb, double beta, double* C, Int ldc,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dgemm(CblasColMajor, TRANS(transa), TRANS(transb), m, n, k,
                           alpha, A, lda, B, ldb, beta, C, ldc);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_gemm);
}

void callAndTime_dsyrk(char uplo, char trans, Int n, Int k, double alpha,
                       const double* A, Int lda, double beta, double* C,
                       Int ldc, DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dsyrk(CblasColMajor, UPLO(uplo), TRANS(trans), n, k, alpha,
                           A, lda, beta, C, ldc);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_syrk);
}

void callAndTime_dtrsm(char side, char uplo, char trans, char diag, Int m,
                       Int n, double alpha, const double* A, Int lda, double* B,
                       Int ldb, DataCollector& data) {
  HIPO_CLOCK_CREATE;
  HighsExtras::blas::dtrsm(CblasColMajor, SIDE(side), UPLO(uplo), TRANS(trans),
                           DIAG(diag), m, n, alpha, A, lda, B, ldb);
  HIPO_CLOCK_STOP(3, data, kTimeBlas_trsm);
}

}  // namespace hipo
