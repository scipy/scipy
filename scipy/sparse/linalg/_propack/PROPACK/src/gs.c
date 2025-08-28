#include "gs.h"


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
void smgs(int n, int k, float* V, int ldv, float* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    /**
     * PROPACK encodes sentinels in 1-index specific way.
     * Therefore, we have to guard for a few edge cases of 0-indexing.
     * TODO: Fix this by rewriting the indexing logic.
     */
    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-based indexing)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)' * vnew
            float coef = sdot_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            float neg_coef = -coef;
            saxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void dmgs(int n, int k, double* V, int ldv, double* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    /**
     * PROPACK encodes sentinels in 1-index specific way.
     * Therefore, we have to guard for a few edge cases of 0-indexing.
     * TODO: Fix this by rewriting the indexing logic.
     */
    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-based indexing)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)' * vnew
            double coef = ddot_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            double neg_coef = -coef;
            daxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void cmgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-based indexing)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)^H * vnew (conjugate dot product)
            PROPACK_CPLXF_TYPE coef = cdotc_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            PROPACK_CPLXF_TYPE neg_coef = PROPACK_cplxf(-crealf(coef), -cimagf(coef));
            caxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void zmgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-based indexing)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)^H * vnew (conjugate dot product)
            PROPACK_CPLX_TYPE coef = zdotc_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            PROPACK_CPLX_TYPE neg_coef = PROPACK_cplx(-creal(coef), -cimag(coef));
            zaxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void scgs(int n, int k, float* V, int ldv, float* vnew, const int* indices, float* work) {
    int ione = 1;
    float one = 1.0f;
    float zero = 0.0f;
    float neg_one = -1.0f;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^T * vnew
        sgemv_("T", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        sgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void dcgs(int n, int k, double* V, int ldv, double* vnew, const int* indices, double* work) {
    int ione = 1;
    double one = 1.0;
    double zero = 0.0;
    double neg_one = -1.0;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^T * vnew
        dgemv_("T", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        dgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void ccgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices, PROPACK_CPLXF_TYPE* work) {
    int ione = 1;
    PROPACK_CPLXF_TYPE one = PROPACK_cplxf(1.0f, 0.0f);
    PROPACK_CPLXF_TYPE zero = PROPACK_cplxf(0.0f, 0.0f);
    PROPACK_CPLXF_TYPE neg_one = PROPACK_cplxf(-1.0f, 0.0f);

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^H * vnew
        cgemv_("C", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        cgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void zcgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices, PROPACK_CPLX_TYPE* work) {
    int ione = 1;
    PROPACK_CPLX_TYPE one = PROPACK_cplx(1.0, 0.0);
    PROPACK_CPLX_TYPE zero = PROPACK_cplx(0.0, 0.0);
    PROPACK_CPLX_TYPE neg_one = PROPACK_cplx(-1.0, 0.0);

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^H * vnew
        zgemv_("C", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        zgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


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
void sreorth(int n, int k, float* V, int ldv, float* vnew, float* normvnew, const int* indices, float alpha, float* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (int itry = 0; itry < NTRY; itry++)
    {
        float normvnew_0 = *normvnew;

        if (iflag == 1) {
            scgs(n, k, V, ldv, vnew, indices, work);
        } else {
            smgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = snrm2_(&n, vnew, &ione);
        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0f;
    for (int i = 0; i < n; i++) { vnew[i] = 0.0f; }
}


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
void dreorth(int n, int k, double* V, int ldv, double* vnew, double* normvnew, const int* indices, double alpha, double* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }
    for (int itry = 0; itry < NTRY; itry++)
    {
        double normvnew_0 = *normvnew;

        if (iflag == 1) {
            dcgs(n, k, V, ldv, vnew, indices, work);
        } else {
            dmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = dnrm2_(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0;
    for (int i = 0; i < n; i++) { vnew[i] = 0.0; }
}


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
void creorth(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, float* normvnew, const int* indices, float alpha, PROPACK_CPLXF_TYPE* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (int itry = 0; itry < NTRY; itry++)
    {
        float normvnew_0 = *normvnew;

        if (iflag == 1) {
            ccgs(n, k, V, ldv, vnew, indices, work);
        } else {
            cmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = scnrm2_(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0f;
    for (int i = 0; i < n; i++) { vnew[i] = PROPACK_cplxf(0.0f, 0.0f); }
}


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
void zreorth(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, double* normvnew, const int* indices, double alpha, PROPACK_CPLX_TYPE* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (int itry = 0; itry < NTRY; itry++)
    {
        double normvnew_0 = *normvnew;

        if (iflag == 1) {
            zcgs(n, k, V, ldv, vnew, indices, work);
        } else {
            zmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = dznrm2_(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0;
    for (int i = 0; i < n; i++) { vnew[i] = PROPACK_cplx(0.0, 0.0); }
}
