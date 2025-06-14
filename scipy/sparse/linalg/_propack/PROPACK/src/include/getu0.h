#ifndef PROPACK_GETU0_H
#define PROPACK_GETU0_H

#include <math.h>
#include <complex.h>
#include "types.h"
#include "blaslapack_declarations.h"
#include "gs.h"
#include "common.h"

/**
 * Generate random vector in span(Op(A)) orthogonal to span(U) - Single precision
 *
 * Attempts to generate a pseudo-random vector in SPAN(Op(A)) orthogonal to
 * span(U(:,0:j-1)), where Op(A) = A if transa=0 and Op(A) = A^T if transa=1.
 *
 * @param transa     0 for A*x, 1 for A^T*x
 * @param m          Number of rows in A
 * @param n          Number of columns in A
 * @param j          Number of existing vectors in U to orthogonalize against
 * @param ntry       Maximum number of attempts
 * @param u0         Output vector (modified in-place)
 * @param u0norm     Pointer to norm of u0 (output)
 * @param U          Matrix of existing vectors (column-major, ldu x j)
 * @param ldu        Leading dimension of U
 * @param aprod      Matrix-vector operation callback function
 * @param dparm      Double parameters for aprod
 * @param iparm      Integer parameters for aprod
 * @param ierr       Error flag (0=success, -1=failure)
 * @param icgs       Orthogonalization method (0=MGS, 1=CGS)
 * @param anormest   Estimate of operator norm (output)
 * @param work       Work array of size max(m,n)
 * @param rng_state  User-supplied random number generator state (uint64_t[4])
 */
void sgetu0(int transa, int m, int n, int j, int ntry, float* u0, float* u0norm, float* U, int ldu,
           PROPACK_aprod_s aprod, float* dparm, int* iparm, int* ierr, int icgs, float* anormest, float* work, uint64_t* rng_state);

/**
 * Generate random vector in span(Op(A)) orthogonal to span(U) - Double precision
 *
 * Attempts to generate a pseudo-random vector in SPAN(Op(A)) orthogonal to
 * span(U(:,0:j-1)), where Op(A) = A if transa=0 and Op(A) = A^T if transa=1.
 *
 * @param transa     0 for A*x, 1 for A^T*x
 * @param m          Number of rows in A
 * @param n          Number of columns in A
 * @param j          Number of existing vectors in U to orthogonalize against
 * @param ntry       Maximum number of attempts
 * @param u0         Output vector (modified in-place)
 * @param u0norm     Pointer to norm of u0 (output)
 * @param U          Matrix of existing vectors (column-major, ldu x j)
 * @param ldu        Leading dimension of U
 * @param aprod      Matrix-vector operation callback function
 * @param dparm      Double parameters for aprod
 * @param iparm      Integer parameters for aprod
 * @param ierr       Error flag (0=success, -1=failure)
 * @param icgs       Orthogonalization method (0=MGS, 1=CGS)
 * @param anormest   Estimate of operator norm (output)
 * @param work       Work array of size max(m,n)
 * @param rng_state  User-supplied random number generator state (uint64_t[4])
 */
void dgetu0(int transa, int m, int n, int j, int ntry, double* u0, double* u0norm, double* U, int ldu,
           PROPACK_aprod_d aprod, double* dparm, int* iparm, int* ierr, int icgs, double* anormest, double* work, uint64_t* rng_state);

/**
 * Generate random vector in span(Op(A)) orthogonal to span(U) - Single precision complex
 *
 * Attempts to generate a pseudo-random vector in SPAN(Op(A)) orthogonal to
 * span(U(:,0:j-1)), where Op(A) = A if transa=0 and Op(A) = A^H if transa=1.
 *
 * @param transa     0 for A*x, 1 for A^H*x
 * @param m          Number of rows in A
 * @param n          Number of columns in A
 * @param j          Number of existing vectors in U to orthogonalize against
 * @param ntry       Maximum number of attempts
 * @param u0         Output vector (modified in-place)
 * @param u0norm     Pointer to norm of u0 (output)
 * @param U          Matrix of existing vectors (column-major, ldu x j)
 * @param ldu        Leading dimension of U
 * @param aprod      Matrix-vector operation callback function
 * @param cparm      Float complex parameters for aprod
 * @param iparm      Integer parameters for aprod
 * @param ierr       Error flag (0=success, -1=failure)
 * @param icgs       Orthogonalization method (0=MGS, 1=CGS)
 * @param anormest   Estimate of operator norm (output)
 * @param work       Work array of size max(m,n)
 * @param rng_state  User-supplied random number generator state (uint64_t[4])
 */
void cgetu0(int transa, int m, int n, int j, int ntry, PROPACK_CPLXF_TYPE* u0, float* u0norm, PROPACK_CPLXF_TYPE* U, int ldu,
           PROPACK_aprod_c aprod, PROPACK_CPLXF_TYPE* cparm, int* iparm, int* ierr, int icgs, float* anormest, PROPACK_CPLXF_TYPE* work,
           uint64_t* rng_state);

/**
 * Generate random vector in span(Op(A)) orthogonal to span(U) - Double precision complex
 *
 * Attempts to generate a pseudo-random vector in SPAN(Op(A)) orthogonal to
 * span(U(:,0:j-1)), where Op(A) = A if transa=0 and Op(A) = A^H if transa=1.
 *
 * @param transa     0 for A*x, 1 for A^H*x
 * @param m          Number of rows in A
 * @param n          Number of columns in A
 * @param j          Number of existing vectors in U to orthogonalize against
 * @param ntry       Maximum number of attempts
 * @param u0         Output vector (modified in-place)
 * @param u0norm     Pointer to norm of u0 (output)
 * @param U          Matrix of existing vectors (column-major, ldu x j)
 * @param ldu        Leading dimension of U
 * @param aprod      Matrix-vector operation callback function
 * @param zparm      Double complex parameters for aprod
 * @param iparm      Integer parameters for aprod
 * @param ierr       Error flag (0=success, -1=failure)
 * @param icgs       Orthogonalization method (0=MGS, 1=CGS)
 * @param anormest   Estimate of operator norm (output)
 * @param work       Work array of size max(m,n)
 * @param rng_state  User-supplied random number generator state (uint64_t[4])
 */
void zgetu0(int transa, int m, int n, int j, int ntry, PROPACK_CPLX_TYPE* u0, double* u0norm, PROPACK_CPLX_TYPE* U, int ldu,
           PROPACK_aprod_z aprod, PROPACK_CPLX_TYPE* zparm, int* iparm, int* ierr, int icgs, double* anormest, PROPACK_CPLX_TYPE* work,
           uint64_t* rng_state);


#endif
