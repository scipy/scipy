#ifndef PROPACK_GEMM_OVWR_H
#define PROPACK_GEMM_OVWR_H

#include "blaslapack_declarations.h"


/**
 * Compute B <- alpha*op(A)*B + beta*B using blocked algorithm
 *
 * @param transa    0 for A, 1 for A^T
 * @param m         number of rows in op(A) and B
 * @param n         number of columns in B
 * @param k         number of columns in op(A)
 * @param alpha     scalar multiplier for A*B
 * @param A         input matrix A (m x k if transa='N', k x m if transa='T')
 * @param lda       leading dimension of A
 * @param beta      scalar multiplier for B
 * @param B         input/output matrix B (m x n), overwritten with result
 * @param ldb       leading dimension of B
 * @param work      workspace array of size at least m * blocksize
 * @param blocksize number of columns to process at once
 */
void sgemm_ovwr(const int transa, int m, int n, int k, float alpha, float* restrict A, int lda, float beta, float* restrict B, int ldb, float* restrict work, int blocksize);


/**
 * Compute A <- alpha*A*op(B) using blocked algorithm
 *
 * @param transb    0 for B, 1 for B^T
 * @param m         number of rows in A
 * @param n         number of columns in A and rows in op(B)
 * @param k         number of columns in op(B)
 * @param alpha     scalar multiplier for A*op(B)
 * @param A         input/output matrix A (m x n), overwritten with result
 * @param lda       leading dimension of A
 * @param B         input matrix B (n x k if transb='N', k x n if transb='T')
 * @param ldb       leading dimension of B
 * @param work      workspace array of size at least blocksize * n
 * @param blocksize number of rows to process at once
 */
void sgemm_ovwr_left(const int transb, int m, int n, int k, float alpha, float* restrict A, int lda, float* restrict B, int ldb, float* restrict work, int blocksize);


/**
 * Compute B <- alpha*op(A)*B + beta*B using blocked algorithm
 *
 * @param transa    0 for A, 1 for A^T
 * @param m         number of rows in op(A) and B
 * @param n         number of columns in B
 * @param k         number of columns in op(A)
 * @param alpha     scalar multiplier for A*B
 * @param A         input matrix A (m x k if transa='N', k x m if transa='T')
 * @param lda       leading dimension of A
 * @param beta      scalar multiplier for B
 * @param B         input/output matrix B (m x n), overwritten with result
 * @param ldb       leading dimension of B
 * @param work      workspace array of size at least m * blocksize
 * @param blocksize number of columns to process at once
 */
void dgemm_ovwr(const int transa, int m, int n, int k, double alpha, double* restrict A, int lda, double beta, double* restrict B, int ldb, double* restrict work, int blocksize);


/**
 * Compute A <- alpha*A*op(B) using blocked algorithm
 *
 * @param transb    0 for B, 1 for B^T
 * @param m         number of rows in A
 * @param n         number of columns in A and rows in op(B)
 * @param k         number of columns in op(B)
 * @param alpha     scalar multiplier for A*op(B)
 * @param A         input/output matrix A (m x n), overwritten with result
 * @param lda       leading dimension of A
 * @param B         input matrix B (n x k if transb='N', k x n if transb='T')
 * @param ldb       leading dimension of B
 * @param work      workspace array of size at least blocksize * n
 * @param blocksize number of rows to process at once
 */
void dgemm_ovwr_left(const int transb, int m, int n, int k, double alpha, double* restrict A, int lda, double* restrict B, int ldb, double* restrict work, int blocksize);


/**
 * Compute A <- A*op(B) using blocked algorithm (A complex, B real)
 *
 * @param transb    0 for B, 1 for B^T
 * @param m         number of rows in A
 * @param n         number of columns in A and rows in op(B)
 * @param k         number of columns in op(B)
 * @param A         input/output complex matrix A (m x n), overwritten with result
 * @param lda       leading dimension of A
 * @param B         input real matrix B (n x k if transb=0, k x n if transb=1)
 * @param ldb       leading dimension of B
 * @param work      workspace array of size at least blocksize * n
 * @param blocksize number of rows to process at once
 */
void csgemm_ovwr_left(const int transb, int m, int n, int k, complex float* restrict A, int lda, const float* restrict B, int ldb, complex float* restrict work, int blocksize);


/**
 * Compute A <- A*op(B) using blocked algorithm (A complex, B real)
 *
 * @param transb    0 for B, 1 for B^T
 * @param m         number of rows in A
 * @param n         number of columns in A and rows in op(B)
 * @param k         number of columns in op(B)
 * @param A         input/output complex matrix A (m x n), overwritten with result
 * @param lda       leading dimension of A
 * @param B         input real matrix B (n x k if transb=0, k x n if transb=1)
 * @param ldb       leading dimension of B
 * @param work      workspace array of size at least blocksize * n
 * @param blocksize number of rows to process at once
 */
void zdgemm_ovwr_left(const int transb, int m, int n, int k, complex double* restrict A, int lda, const double* restrict B, int ldb, complex double* restrict work, int blocksize);

#endif
