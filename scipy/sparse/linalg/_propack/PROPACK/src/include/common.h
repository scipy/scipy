#ifndef PROPACK__COMMON_H
#define PROPACK__COMMON_H
#include <stdint.h>
#include "blaslapack_declarations.h"

/**
 * Generate a random single-precision floating point number in [0, 1).
 *
 * @param state Pointer to 4-element uint64_t array containing xoshiro256+ PRNG state
 * @return Random float in the range [-1.0f, 1.0f)
 */
float random_float(uint64_t* state);

/**
 * Generate a random double-precision floating point number in [0, 1).
 *
 * @param state Pointer to 4-element uint64_t array containing xoshiro256+ PRNG state
 * @return Random double in the range [-1.0, 1.0)
 */
double random_double(uint64_t* state);

/**
 * Perform an implicit LQ SVD sweep with shift sigma (single precision).
 *
 * Applies Givens rotations to chase a bulge down a lower bidiagonal matrix,
 * performing one step of the implicit QR algorithm for SVD computation.
 *
 * @param jobu If nonzero, update left singular vectors U
 * @param jobv If nonzero, update right singular vectors V
 * @param m Number of rows in U matrix
 * @param n Number of rows in V matrix
 * @param k Size of the bidiagonal submatrix to process
 * @param sigma Shift value for the implicit QR step
 * @param D Diagonal elements of the bidiagonal matrix (length k)
 * @param E Superdiagonal elements of the bidiagonal matrix (length k-1)
 * @param U Left singular vectors matrix (m x ldu), updated if jobu != 0
 * @param ldu Leading dimension of U
 * @param V Right singular vectors matrix (n x ldv), updated if jobv != 0
 * @param ldv Leading dimension of V
 */
void sbsvdstep(const int jobu, const int jobv, int m, int n, int k, float sigma, float* D, float* E, float* U, int ldu, float* V, int ldv);

/**
 * Perform an implicit LQ SVD sweep with shift sigma (double precision).
 *
 * Applies Givens rotations to chase a bulge down a lower bidiagonal matrix,
 * performing one step of the implicit QR algorithm for SVD computation.
 *
 * @param jobu If nonzero, update left singular vectors U
 * @param jobv If nonzero, update right singular vectors V
 * @param m Number of rows in U matrix
 * @param n Number of rows in V matrix
 * @param k Size of the bidiagonal submatrix to process
 * @param sigma Shift value for the implicit QR step
 * @param D Diagonal elements of the bidiagonal matrix (length k)
 * @param E Superdiagonal elements of the bidiagonal matrix (length k-1)
 * @param U Left singular vectors matrix (m x ldu), updated if jobu != 0
 * @param ldu Leading dimension of U
 * @param V Right singular vectors matrix (n x ldv), updated if jobv != 0
 * @param ldv Leading dimension of V
 */
void dbsvdstep(const int jobu, const int jobv, int m, int n, int k, double sigma, double* D, double* E, double* U, int ldu, double* V, int ldv);

/**
 * QR factorization of a lower bidiagonal matrix using Givens rotations (single precision).
 *
 * Computes the QR factorization of a lower bidiagonal matrix, optionally
 * accumulating the orthogonal transformation matrix Qt.
 *
 * @param ignorelast If nonzero, skip processing the last element
 * @param jobq If nonzero, accumulate orthogonal transformations in Qt
 * @param n Size of the bidiagonal matrix
 * @param D Diagonal elements (length n), modified in-place
 * @param E Subdiagonal elements (length n-1), modified in-place
 * @param c1 Output: cosine of final Givens rotation
 * @param c2 Output: sine of final Givens rotation
 * @param Qt Orthogonal transformation matrix (n+1 x ldq), updated if jobq != 0
 * @param ldq Leading dimension of Qt
 */
void sbdqr(const int ignorelast, const int jobq, const int n, float* restrict D, float* restrict E, float* c1, float* c2, float* restrict Qt, int ldq);

/**
 * QR factorization of a lower bidiagonal matrix using Givens rotations (double precision).
 *
 * Computes the QR factorization of a lower bidiagonal matrix, optionally
 * accumulating the orthogonal transformation matrix Qt.
 *
 * @param ignorelast If nonzero, skip processing the last element
 * @param jobq If nonzero, accumulate orthogonal transformations in Qt
 * @param n Size of the bidiagonal matrix
 * @param D Diagonal elements (length n), modified in-place
 * @param E Subdiagonal elements (length n-1), modified in-place
 * @param c1 Output: cosine of final Givens rotation
 * @param c2 Output: sine of final Givens rotation
 * @param Qt Orthogonal transformation matrix (n+1 x ldq), updated if jobq != 0
 * @param ldq Leading dimension of Qt
 */
void dbdqr(const int ignorelast, const int jobq, const int n, double* restrict D, double* restrict E, double* c1, double* c2, double* restrict Qt, int ldq);

/**
 * Refine error bounds on Ritz values by clustering and gap analysis (single precision).
 *
 * Improves error bounds by identifying clusters of close Ritz values and
 * computing gaps between eigenvalues for better convergence estimates.
 *
 * @param n Total number of eigenvalues
 * @param k Number of computed Ritz values
 * @param theta Array of Ritz values (length k)
 * @param bound Error bounds for Ritz values (length k), modified in-place
 * @param tol Tolerance for bound refinement
 * @param eps34 Clustering threshold (typically machine precision^(3/4))
 */
void srefinebounds(const int n, const int k, float* restrict theta, float* restrict bound, const float tol, const float eps34);

/**
 * Refine error bounds on Ritz values by clustering and gap analysis (double precision).
 *
 * Improves error bounds by identifying clusters of close Ritz values and
 * computing gaps between eigenvalues for better convergence estimates.
 *
 * @param n Total number of eigenvalues
 * @param k Number of computed Ritz values
 * @param theta Array of Ritz values (length k)
 * @param bound Error bounds for Ritz values (length k), modified in-place
 * @param tol Tolerance for bound refinement
 * @param eps34 Clustering threshold (typically machine precision^(3/4))
 */
void drefinebounds(const int n, const int k, double* restrict theta, double* restrict bound, const double tol, const double eps34);

/**
 * Find intervals where |mu| exceeds thresholds for reorthogonalization (single precision).
 *
 * Identifies peaks in the mu array that exceed delta and finds the intervals
 * around these peaks where values remain above eta, used for selective
 * reorthogonalization in Lanczos-type algorithms.
 *
 * @param mu Array of reorthogonalization indicators (length j)
 * @param j Length of the mu array
 * @param delta Peak detection threshold
 * @param eta Interval boundary threshold
 * @param indices Output array of interval boundaries, terminated by j
 */
void scompute_int(float* restrict mu, const int j, const float delta, const float eta, int* restrict indices);

/**
 * Set mu values to a constant within specified intervals (single precision).
 *
 * @param k Index up to which to process intervals
 * @param mu Array of values to modify
 * @param indices Array of interval boundaries (pairs of start, end indices)
 * @param val Value to assign within the intervals
 */
void sset_mu(const int k, float* restrict mu, int* const restrict indices, const float val);

/**
 * Update mu recurrence relation for Lanczos reorthogonalization (single precision).
 *
 * Computes the recurrence relation for tracking orthogonality loss in
 * Lanczos-type algorithms, updating the mu array and tracking maximum values.
 *
 * @param mumax Output: maximum absolute value in updated mu array
 * @param mu Reorthogonalization tracking array, modified in-place
 * @param nu Auxiliary tracking array
 * @param j Current iteration index
 * @param alpha Diagonal elements of tridiagonal matrix
 * @param beta Off-diagonal elements of tridiagonal matrix
 * @param anorm Norm of the input matrix
 * @param eps1 Numerical precision parameter
 */
void supdate_mu(float* mumax, float* restrict mu, float* restrict nu, const int j, float* restrict alpha, float* restrict beta, const float anorm, const float eps1);

/**
 * Update nu recurrence relation for Lanczos reorthogonalization (single precision).
 *
 * Computes the recurrence relation for the nu auxiliary array used in
 * tracking orthogonality loss in Lanczos-type algorithms.
 *
 * @param numax Output: maximum absolute value in updated nu array
 * @param mu Reorthogonalization tracking array
 * @param nu Auxiliary tracking array, modified in-place
 * @param j Current iteration index
 * @param alpha Diagonal elements of tridiagonal matrix
 * @param beta Off-diagonal elements of tridiagonal matrix
 * @param anorm Norm of the input matrix
 * @param eps1 Numerical precision parameter
 */
void supdate_nu(float* numax, float* restrict mu, float* restrict nu, const int j, float* restrict alpha, float* restrict beta, const float anorm, const float eps1);

/**
 * Find intervals where |mu| exceeds thresholds for reorthogonalization (double precision).
 *
 * Identifies peaks in the mu array that exceed delta and finds the intervals
 * around these peaks where values remain above eta, used for selective
 * reorthogonalization in Lanczos-type algorithms.
 *
 * @param mu Array of reorthogonalization indicators (length j)
 * @param j Length of the mu array
 * @param delta Peak detection threshold
 * @param eta Interval boundary threshold
 * @param indices Output array of interval boundaries, terminated by j
 */
void dcompute_int(double* restrict mu, const int j, const double delta, const double eta, int* restrict indices);

/**
 * Set mu values to a constant within specified intervals (double precision).
 *
 * @param k Index up to which to process intervals
 * @param mu Array of values to modify
 * @param indices Array of interval boundaries (pairs of start, end indices)
 * @param val Value to assign within the intervals
 */
void dset_mu(const int k, double* restrict mu, int* const restrict indices, const double val);

/**
 * Update mu recurrence relation for Lanczos reorthogonalization (double precision).
 *
 * Computes the recurrence relation for tracking orthogonality loss in
 * Lanczos-type algorithms, updating the mu array and tracking maximum values.
 *
 * @param mumax Output: maximum absolute value in updated mu array
 * @param mu Reorthogonalization tracking array, modified in-place
 * @param nu Auxiliary tracking array
 * @param j Current iteration index
 * @param alpha Diagonal elements of tridiagonal matrix
 * @param beta Off-diagonal elements of tridiagonal matrix
 * @param anorm Norm of the input matrix
 * @param eps1 Numerical precision parameter
 */
void dupdate_mu(double* mumax, double* restrict mu, double* restrict nu, const int j, double* restrict alpha, double* restrict beta, const double anorm, const double eps1);

/**
 * Update nu recurrence relation for Lanczos reorthogonalization (double precision).
 *
 * Computes the recurrence relation for the nu auxiliary array used in
 * tracking orthogonality loss in Lanczos-type algorithms.
 *
 * @param numax Output: maximum absolute value in updated nu array
 * @param mu Reorthogonalization tracking array
 * @param nu Auxiliary tracking array, modified in-place
 * @param j Current iteration index
 * @param alpha Diagonal elements of tridiagonal matrix
 * @param beta Off-diagonal elements of tridiagonal matrix
 * @param anorm Norm of the input matrix
 * @param eps1 Numerical precision parameter
 */
void dupdate_nu(double* numax, double* restrict mu, double* restrict nu, const int j, double* restrict alpha, double* restrict beta, const double anorm, const double eps1);


#endif // PROPACK__COMMON_H
