#ifndef PROPACK_GS_H
#define PROPACK_GS_H

#include "blaslapack_declarations.h"


/**
 * @brief Performs Modified Gram-Schmidt (MGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against the `k` columns of the matrix
 * `V` using the Modified Gram-Schmidt process. The orthogonalization is performed
 * iteratively over blocks of columns specified by the `indices` array.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 */
void smgs(int n, int k, float* V, int ldv, float* vnew, const int* indices);


/**
 * @brief Performs Modified Gram-Schmidt (MGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against the `k` columns of the matrix
 * `V` using the Modified Gram-Schmidt process. The orthogonalization is performed
 * iteratively over blocks of columns specified by the `indices` array.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 */
void dmgs(int n, int k, double* V, int ldv, double* vnew, const int* indices);


/**
 * @brief Performs Modified Gram-Schmidt (MGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against the `k` columns of the matrix
 * `V` using the Modified Gram-Schmidt process. The orthogonalization is performed
 * iteratively over blocks of columns specified by the `indices` array.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 */
void cmgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices);


/**
 * @brief Performs Modified Gram-Schmidt (MGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against the `k` columns of the matrix
 * `V` using the Modified Gram-Schmidt process. The orthogonalization is performed
 * iteratively over blocks of columns specified by the `indices` array.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 */
void zmgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices);


/**
 * @brief Performs block Classical Gram-Schmidt (CGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against blocks of columns in the matrix `V`
 * using the Classical Gram-Schmidt process. The blocks of columns are specified by the
 * `indices` array, and the orthogonalization is performed iteratively for each block.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 * @param work    Pointer to a workspace array of length at least equal to the size of the
 *                largest block (end - start + 1) in the `indices` array.
 */
void scgs(int n, int k, float* V, int ldv, float* vnew, const int* indices, float* work);


/**
 * @brief Performs block Classical Gram-Schmidt (CGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against blocks of columns in the matrix `V`
 * using the Classical Gram-Schmidt process. The blocks of columns are specified by the
 * `indices` array, and the orthogonalization is performed iteratively for each block.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 * @param work    Pointer to a workspace array of length at least equal to the size of the
 *                largest block (end - start + 1) in the `indices` array.
 */
void dcgs(int n, int k, double* V, int ldv, double* vnew, const int* indices, double* work);


/**
 * @brief Performs block Classical Gram-Schmidt (CGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against blocks of columns in the matrix `V`
 * using the Classical Gram-Schmidt process. The blocks of columns are specified by the
 * `indices` array, and the orthogonalization is performed iteratively for each block.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 * @param work    Pointer to a workspace array of length at least equal to the size of the
 *                largest block (end - start + 1) in the `indices` array.
 */
void ccgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices, PROPACK_CPLXF_TYPE* work);


/**
 * @brief Performs block Classical Gram-Schmidt (CGS) orthogonalization on a vector.
 *
 * This function orthogonalizes the vector `vnew` against blocks of columns in the matrix `V`
 * using the Classical Gram-Schmidt process. The blocks of columns are specified by the
 * `indices` array, and the orthogonalization is performed iteratively for each block.
 *
 * @param n       The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k       The number of columns in the matrix `V` to orthogonalize against.
 * @param V       Pointer to the matrix `V` stored in column-major order.
 * @param ldv     The leading dimension of the matrix `V`.
 * @param vnew    Pointer to the vector to be orthogonalized. Modified in-place.
 * @param indices Pointer to an array specifying the start and end indices of column blocks
 *                in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                for each block. The last element is marked as a terminator typically a value
 *                greater than k-1.
 * @param work    Pointer to a workspace array of length at least equal to the size of the
 *                largest block (end - start + 1) in the `indices` array.
 */
void zcgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices, PROPACK_CPLX_TYPE* work);


/**
 * @brief Reorthogonalizes a vector against a subset of columns of a matrix.
 *
 * This function reorthogonalizes the vector `vnew` against a subset of the columns
 * of the matrix `V` using iterated Classical or Modified Gram-Schmidt. The process
 * is repeated until the norm of the reorthogonalized vector satisfies the condition:
 *
 *     ||vnew'|| > alpha * ||vnew||
 *
 * If the condition is not satisfied after a fixed number of attempts, the vector
 * `vnew` is deemed to lie numerically in the span of the selected columns of `V`
 * and is set to the zero vector.
 *
 * @param n         The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k         The number of columns in the matrix `V` to orthogonalize against.
 * @param V         Pointer to the matrix `V` stored in column-major order.
 * @param ldv       The leading dimension of the matrix `V`.
 * @param vnew      Pointer to the vector to be reorthogonalized. Modified in-place.
 * @param normvnew  Pointer to the norm of the vector `vnew`. Updated in-place.
 * @param indices   Pointer to an array specifying the start and end indices of column blocks
 *                  in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                  where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                  for each block. The last element is marked as a terminator typically a value
 *                  greater than k-1.
 * @param alpha     The threshold factor for the reorthogonalization condition.
 * @param work      Pointer to a workspace array of length at least equal to the size of the
 *                  largest block (end - start + 1) in the `indices` array. Used only if `iflag == 1`.
 * @param iflag     Determines the orthogonalization method:
 *                  - `0`: Iterated Modified Gram-Schmidt (MGS).
 *                  - `1`: Iterated Classical Gram-Schmidt (CGS).
 */
void sreorth(int n, int k, float* V, int ldv, float* vnew, float* normvnew, const int* indices, float alpha, float* work, int iflag);


/**
 * @brief Reorthogonalizes a vector against a subset of columns of a matrix.
 *
 * This function reorthogonalizes the vector `vnew` against a subset of the columns
 * of the matrix `V` using iterated Classical or Modified Gram-Schmidt. The process
 * is repeated until the norm of the reorthogonalized vector satisfies the condition:
 *
 *     ||vnew'|| > alpha * ||vnew||
 *
 * If the condition is not satisfied after a fixed number of attempts, the vector
 * `vnew` is deemed to lie numerically in the span of the selected columns of `V`
 * and is set to the zero vector.
 *
 * @param n         The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k         The number of columns in the matrix `V` to orthogonalize against.
 * @param V         Pointer to the matrix `V` stored in column-major order.
 * @param ldv       The leading dimension of the matrix `V`.
 * @param vnew      Pointer to the vector to be reorthogonalized. Modified in-place.
 * @param normvnew  Pointer to the norm of the vector `vnew`. Updated in-place.
 * @param indices   Pointer to an array specifying the start and end indices of column blocks
 *                  in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                  where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                  for each block. The last element is marked as a terminator typically a value
 *                  greater than k-1.
 * @param alpha     The threshold factor for the reorthogonalization condition.
 * @param work      Pointer to a workspace array of length at least equal to the size of the
 *                  largest block (end - start + 1) in the `indices` array. Used only if `iflag == 1`.
 * @param iflag     Determines the orthogonalization method:
 *                  - `0`: Iterated Modified Gram-Schmidt (MGS).
 *                  - `1`: Iterated Classical Gram-Schmidt (CGS).
 */
void dreorth(int n, int k, double* V, int ldv, double* vnew, double* normvnew, const int* indices, double alpha, double* work, int iflag);


/**
 * @brief Reorthogonalizes a vector against a subset of columns of a matrix.
 *
 * This function reorthogonalizes the vector `vnew` against a subset of the columns
 * of the matrix `V` using iterated Classical or Modified Gram-Schmidt. The process
 * is repeated until the norm of the reorthogonalized vector satisfies the condition:
 *
 *     ||vnew'|| > alpha * ||vnew||
 *
 * If the condition is not satisfied after a fixed number of attempts, the vector
 * `vnew` is deemed to lie numerically in the span of the selected columns of `V`
 * and is set to the zero vector.
 *
 * @param n         The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k         The number of columns in the matrix `V` to orthogonalize against.
 * @param V         Pointer to the matrix `V` stored in column-major order.
 * @param ldv       The leading dimension of the matrix `V`.
 * @param vnew      Pointer to the vector to be reorthogonalized. Modified in-place.
 * @param normvnew  Pointer to the norm of the vector `vnew`. Updated in-place.
 * @param indices   Pointer to an array specifying the start and end indices of column blocks
 *                  in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                  where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                  for each block. The last element is marked as a terminator typically a value
 *                  greater than k-1.
 * @param alpha     The threshold factor for the reorthogonalization condition.
 * @param work      Pointer to a workspace array of length at least equal to the size of the
 *                  largest block (end - start + 1) in the `indices` array. Used only if `iflag == 1`.
 * @param iflag     Determines the orthogonalization method:
 *                  - `0`: Iterated Modified Gram-Schmidt (MGS).
 *                  - `1`: Iterated Classical Gram-Schmidt (CGS).
 */
void creorth(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, float* normvnew, const int* indices, float alpha, PROPACK_CPLXF_TYPE* work, int iflag);


/**
 * @brief Reorthogonalizes a vector against a subset of columns of a matrix.
 *
 * This function reorthogonalizes the vector `vnew` against a subset of the columns
 * of the matrix `V` using iterated Classical or Modified Gram-Schmidt. The process
 * is repeated until the norm of the reorthogonalized vector satisfies the condition:
 *
 *     ||vnew'|| > alpha * ||vnew||
 *
 * If the condition is not satisfied after a fixed number of attempts, the vector
 * `vnew` is deemed to lie numerically in the span of the selected columns of `V`
 * and is set to the zero vector.
 *
 * @param n         The number of rows in the matrix `V` and the length of the vector `vnew`.
 * @param k         The number of columns in the matrix `V` to orthogonalize against.
 * @param V         Pointer to the matrix `V` stored in column-major order.
 * @param ldv       The leading dimension of the matrix `V`.
 * @param vnew      Pointer to the vector to be reorthogonalized. Modified in-place.
 * @param normvnew  Pointer to the norm of the vector `vnew`. Updated in-place.
 * @param indices   Pointer to an array specifying the start and end indices of column blocks
 *                  in `V` to process. The array should have at least `2 * num_blocks` elements,
 *                  where `num_blocks` is the number of blocks. It encodes the pairs [start, end]
 *                  for each block. The last element is marked as a terminator typically a value
 *                  greater than k-1.
 * @param alpha     The threshold factor for the reorthogonalization condition.
 * @param work      Pointer to a workspace array of length at least equal to the size of the
 *                  largest block (end - start + 1) in the `indices` array. Used only if `iflag == 1`.
 * @param iflag     Determines the orthogonalization method:
 *                  - `0`: Iterated Modified Gram-Schmidt (MGS).
 *                  - `1`: Iterated Classical Gram-Schmidt (CGS).
 */
void zreorth(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, double* normvnew, const int* indices, double alpha, PROPACK_CPLX_TYPE* work, int iflag);


#endif /* GS_H */
