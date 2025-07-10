#ifndef PROPACK_GS_H
#define PROPACK_GS_H

#include "blaslapack_declarations.h"


/**
 * Modified Gram-Schmidt orthogonalization - Single precision
 *
 * Orthogonalizes vnew against specified blocks of vectors in V, following
 * the PROPACK partial reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive
 *                 Terminated when start_i <= 0 or start_i > k or start_i > end_i
 */
void smgs(int n, int k, float* V, int ldv, float* vnew, const int* indices);

/**
 * Modified Gram-Schmidt orthogonalization - Double precision
 *
 * Orthogonalizes vnew against specified blocks of vectors in V, following
 * the PROPACK partial reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive
 *                 Terminated when start_i <= 0 or start_i > k or start_i > end_i
 */
void dmgs(int n, int k, double* V, int ldv, double* vnew, const int* indices);

/**
 * Modified Gram-Schmidt orthogonalization - Single precision complex
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using conjugate
 * inner product, following the PROPACK partial reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive
 *                 Terminated when start_i <= 0 or start_i > k or start_i > end_i
 */
void cmgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices);

/**
 * Modified Gram-Schmidt orthogonalization - Double precision complex
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using conjugate
 * inner product, following the PROPACK partial reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive
 *                 Terminated when start_i <= 0 or start_i > k or start_i > end_i
 */
void zmgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices);

/**
 * Classical Gram-Schmidt orthogonalization - Single precision
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * block matrix-vector operations, following the PROPACK partial
 * reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param work     Work array
 */
void scgs(int n, int k, float* V, int ldv, float* vnew, const int* indices, float* work);

/**
 * Classical Gram-Schmidt orthogonalization - Double precision
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * block matrix-vector operations, following the PROPACK partial
 * reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param work     Work array
 */
void dcgs(int n, int k, double* V, int ldv, double* vnew, const int* indices, double* work);

/**
 * Classical Gram-Schmidt orthogonalization - Single precision complex
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * block matrix-vector operations with conjugate transpose, following
 * the PROPACK partial reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param work     Work array
 */
void ccgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices, PROPACK_CPLXF_TYPE* work);

/**
 * Classical Gram-Schmidt orthogonalization - Double precision complex
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * block matrix-vector operations with conjugate transpose, following
 * the PROPACK partial reorthogonalization pattern.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param work     Work array
 */
void zcgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices, PROPACK_CPLX_TYPE* work);

/**
 * Iterated orthogonalization with convergence check - Single precision
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * iterated classical or modified Gram-Schmidt until convergence criterion
 * is met or maximum iterations reached.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param normvnew Pointer to norm of vnew (input: initial norm, output: final norm)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param alpha    Convergence parameter (typically 0.5 or 0.717)
 * @param work     Work array of size at least max(end_i-start_i+1) (only used if iflag==1)
 * @param iflag    0 for Modified Gram-Schmidt, 1 for Classical Gram-Schmidt
 */
void sreorth(int n, int k, float* V, int ldv, float* vnew, float* normvnew, const int* indices, float alpha, float* work, int iflag);

/**
 * Iterated orthogonalization with convergence check - Double precision
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * iterated classical or modified Gram-Schmidt until convergence criterion
 * is met or maximum iterations reached.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param normvnew Pointer to norm of vnew (input: initial norm, output: final norm)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param alpha    Convergence parameter (typically 0.5 or 0.717)
 * @param work     Work array of size at least max(end_i-start_i+1) (only used if iflag==1)
 * @param iflag    0 for Modified Gram-Schmidt, 1 for Classical Gram-Schmidt
 */
void dreorth(int n, int k, double* V, int ldv, double* vnew, double* normvnew, const int* indices, double alpha, double* work, int iflag);

/**
 * Iterated orthogonalization with convergence check - Single precision complex
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * iterated classical or modified Gram-Schmidt until convergence criterion
 * is met or maximum iterations reached.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param normvnew Pointer to norm of vnew (input: initial norm, output: final norm)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param alpha    Convergence parameter (typically 0.5 or 0.717)
 * @param work     Work array of size at least max(end_i-start_i+1) (only used if iflag==1)
 * @param iflag    0 for Modified Gram-Schmidt, 1 for Classical Gram-Schmidt
 */
void creorth(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, float* normvnew, const int* indices, float alpha, PROPACK_CPLXF_TYPE* work, int iflag);

/**
 * Iterated orthogonalization with convergence check - Double precision complex
 *
 * Orthogonalizes vnew against specified blocks of vectors in V using
 * iterated classical or modified Gram-Schmidt until convergence criterion
 * is met or maximum iterations reached.
 *
 * @param n        Length of vectors
 * @param k        Number of vectors in V (maximum column index)
 * @param V        Matrix of orthogonal vectors (column-major, ldv x k)
 * @param ldv      Leading dimension of V
 * @param vnew     Vector to be orthogonalized (modified in-place)
 * @param normvnew Pointer to norm of vnew (input: initial norm, output: final norm)
 * @param indices  Array specifying which columns to orthogonalize against
 *                 Format: [start1, end1, start2, end2, ...]
 *                 Blocks are [start_i, end_i] inclusive (0-based indexing)
 *                 Terminated when start_i < 0 or start_i > k or start_i > end_i
 * @param alpha    Convergence parameter (typically 0.5 or 0.717)
 * @param work     Work array of size at least max(end_i-start_i+1) (only used if iflag==1)
 * @param iflag    0 for Modified Gram-Schmidt, 1 for Classical Gram-Schmidt
 */
void zreorth(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, double* normvnew, const int* indices, double alpha, PROPACK_CPLX_TYPE* work, int iflag);


#endif /* GS_H */
