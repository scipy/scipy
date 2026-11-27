#ifndef FACTORHIGHS_CALL_AND_TIME_BLAS_H
#define FACTORHIGHS_CALL_AND_TIME_BLAS_H

#include "DataCollector.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// level 1
void callAndTime_daxpy(Int n, double da, const double* dx, Int incx, double* dy,
                       Int incy, DataCollector& data);
void callAndTime_dcopy(Int n, const double* dx, Int incx, double* dy, Int incy,
                       DataCollector& data);
void callAndTime_dscal(Int n, const double da, double* dx, Int incx,
                       DataCollector& data);
void callAndTime_dswap(Int n, double* dx, Int incx, double* dy, Int incy,
                       DataCollector& data);

// level 2
void callAndTime_dgemv(char trans, Int m, Int n, double alpha, const double* A,
                       Int lda, const double* x, Int incx, double beta,
                       double* y, Int incy, DataCollector& data);
void callAndTime_dtpsv(char uplo, char trans, char diag, Int n,
                       const double* ap, double* x, Int incx,
                       DataCollector& data);
void callAndTime_dtrsv(char uplo, char trans, char diag, Int n, const double* A,
                       Int lda, double* x, Int incx, DataCollector& data);
void callAndTime_dger(Int m, Int n, double alpha, const double* x, Int incx,
                      const double* y, Int incy, double* A, Int lda,
                      DataCollector& data);

// level 3
void callAndTime_dgemm(char transa, char transb, Int m, Int n, Int k,
                       double alpha, const double* A, Int lda, const double* B,
                       Int ldb, double beta, double* C, Int ldc,
                       DataCollector& data);
void callAndTime_dsyrk(char uplo, char trans, Int n, Int k, double alpha,
                       const double* a, Int lda, double beta, double* c,
                       Int ldc, DataCollector& data);
void callAndTime_dtrsm(char side, char uplo, char trans, char diag, Int m,
                       Int n, double alpha, const double* a, Int lda, double* b,
                       Int ldb, DataCollector& data);

}  // namespace hipo

#endif