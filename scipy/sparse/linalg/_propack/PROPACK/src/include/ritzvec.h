#include "propack/types.h"
#include "common.h"
#include "gemm_overwrite.h"

/**
 * Compute Ritz vectors from bidiagonal singular value decomposition (single precision)
 *
 * Computes singular vectors from the SVD of a bidiagonal matrix obtained from
 * Lanczos bidiagonalization. Uses a two-stage procedure: QR factorization followed
 * by divide-and-conquer SVD.
 *
 * @param which     0 for smallest singular values, 1 for largest singular values
 * @param jobu      1 to compute left singular vectors, 0 otherwise
 * @param jobv      1 to compute right singular vectors, 0 otherwise
 * @param m         number of rows in original matrix
 * @param n         number of columns in original matrix
 * @param k         number of desired singular vectors
 * @param dim       dimension of Krylov subspace
 * @param D         diagonal elements of bidiagonal matrix
 * @param E         super-diagonal elements of bidiagonal matrix
 * @param U         left singular vectors (m x k), modified in-place
 * @param ldu       leading dimension of U
 * @param V         right singular vectors (n x k), modified in-place
 * @param ldv       leading dimension of V
 * @param work      workspace array
 * @param in_lwrk   size of work array
 * @param iwork     integer workspace array
 */
void sritzvec(const int which, const int jobu, const int jobv, const int m, const int n, const int k, int dim,
              float* restrict D, float* restrict E, float* restrict U, const int ldu,
              float* restrict V, const int ldv, float* restrict work, const int in_lwrk, int* restrict iwork);


/**
 * Compute Ritz vectors from bidiagonal singular value decomposition (double precision)
 *
 * Computes singular vectors from the SVD of a bidiagonal matrix obtained from
 * Lanczos bidiagonalization. Uses a two-stage procedure: QR factorization followed
 * by divide-and-conquer SVD.
 *
 * @param which     0 for smallest singular values, 1 for largest singular values
 * @param jobu      1 to compute left singular vectors, 0 otherwise
 * @param jobv      1 to compute right singular vectors, 0 otherwise
 * @param m         number of rows in original matrix
 * @param n         number of columns in original matrix
 * @param k         number of desired singular vectors
 * @param dim       dimension of Krylov subspace
 * @param D         diagonal elements of bidiagonal matrix
 * @param E         super-diagonal elements of bidiagonal matrix
 * @param U         left singular vectors (m x k), modified in-place
 * @param ldu       leading dimension of U
 * @param V         right singular vectors (n x k), modified in-place
 * @param ldv       leading dimension of V
 * @param work      workspace array
 * @param in_lwrk   size of work array
 * @param iwork     integer workspace array
 */
void dritzvec(const int which, const int jobu, const int jobv, const int m, const int n, const int k, int dim,
              double* restrict D, double* restrict E, double* restrict U, const int ldu,
              double* restrict V, const int ldv, double* restrict work, const int in_lwrk, int* restrict iwork);


/**
 * Compute Ritz vectors from bidiagonal SVD (single precision complex)
 *
 * Mixed-precision version where left singular vectors are complex and right singular
 * vectors are real. Computes singular vectors from the SVD of a bidiagonal matrix
 * obtained from Lanczos bidiagonalization.
 *
 * @param which     0 for smallest singular values, 1 for largest singular values
 * @param jobu      1 to compute left singular vectors, 0 otherwise
 * @param jobv      1 to compute right singular vectors, 0 otherwise
 * @param m         number of rows in original matrix
 * @param n         number of columns in original matrix
 * @param k         number of desired singular vectors
 * @param dim       dimension of Krylov subspace
 * @param D         diagonal elements of bidiagonal matrix (real)
 * @param E         super-diagonal elements of bidiagonal matrix (real)
 * @param U         left singular vectors (complex, m x k), modified in-place
 * @param ldu       leading dimension of U
 * @param V         right singular vectors (complex, n x k), modified in-place
 * @param ldv       leading dimension of V
 * @param work      real workspace array
 * @param in_lwrk   size of work array
 * @param cwork     complex workspace array for mixed-precision operations
 * @param lcwrk     size of complex workspace array
 * @param iwork     integer workspace array
 */
void critzvec(const int which, const int jobu, const int jobv, const int m, const int n, const int k, int dim,
              float* restrict D, float* restrict E, PROPACK_CPLXF_TYPE* restrict U, const int ldu,
              PROPACK_CPLXF_TYPE* restrict V, const int ldv, float* restrict work, const int in_lwrk,
              PROPACK_CPLXF_TYPE* restrict cwork, const int lcwrk, int* restrict iwork);


/**
 * Compute Ritz vectors from bidiagonal SVD (double precision complex)
 *
 * Complex version where both left and right singular vectors are complex double.
 * Computes singular vectors from the SVD of a bidiagonal matrix obtained from
 * Lanczos bidiagonalization.
 *
 * @param which     0 for smallest singular values, 1 for largest singular values
 * @param jobu      1 to compute left singular vectors, 0 otherwise
 * @param jobv      1 to compute right singular vectors, 0 otherwise
 * @param m         number of rows in original matrix
 * @param n         number of columns in original matrix
 * @param k         number of desired singular vectors
 * @param dim       dimension of Krylov subspace
 * @param D         diagonal elements of bidiagonal matrix (real)
 * @param E         super-diagonal elements of bidiagonal matrix (real)
 * @param U         left singular vectors (complex, m x k), modified in-place
 * @param ldu       leading dimension of U
 * @param V         right singular vectors (complex, n x k), modified in-place
 * @param ldv       leading dimension of V
 * @param work      real workspace array
 * @param in_lwrk   size of work array
 * @param zwork     complex workspace array for mixed-precision operations
 * @param lzwrk     size of complex workspace array
 * @param iwork     integer workspace array
 */
void zritzvec(const int which, const int jobu, const int jobv, const int m, const int n, const int k, int dim,
              double* restrict D, double* restrict E, PROPACK_CPLX_TYPE* restrict U, const int ldu,
              PROPACK_CPLX_TYPE* restrict V, const int ldv, double* restrict work, const int in_lwrk,
              PROPACK_CPLX_TYPE* restrict zwork, const int lzwrk, int* restrict iwork);
