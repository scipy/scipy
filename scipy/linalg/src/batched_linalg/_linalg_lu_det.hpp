/*
 * Templated loops for `linalg.lu` and `linalg.det`
 */
#pragma once

namespace sp_linalg {

/**
 * @brief Convert a LAPACK 1-based pivot sequence to a 0-based row permutation.
 *
 * On exit, @p perm satisfies P @ A = L @ U: row i of L @ U is row perm[i] of A.
 *
 * @param[in]  ipiv  Pivot sequence of length @p mn from ?getrf, 1-based.
 * @param[out] perm  Permutation of length @p m.
 * @param[in]  m     Number of rows.
 * @param[in]  mn    min(m, n), the length of the pivot sequence.
 */
inline void ipiv_to_perm(const CBLAS_INT *ipiv, CBLAS_INT *perm, const CBLAS_INT m, const CBLAS_INT mn)
{
    // Applying the swaps in reverse order yields the inverse permutation directly.
    for (CBLAS_INT i = 0; i < m; i++) { perm[i] = i; }
    for (CBLAS_INT i = mn - 1; i >= 0; i--) {
        CBLAS_INT tmp = perm[i];
        perm[i] = perm[ipiv[i] - 1];
        perm[ipiv[i] - 1] = tmp;
    }
}


/**
 * @brief Permute the rows of a C-ordered m-by-ncols matrix in place.
 *
 * On exit, row i of @p data is the original row perm[i].
 *
 * @param[in,out] data   Row-major matrix, m * ncols elements.
 * @param[in]     perm   Permutation of length @p m.
 * @param[in]     tmp    Scratch of at least m * ncols elements.
 * @param[in]     m      Number of rows.
 * @param[in]     ncols  Number of columns.
 */
template<typename T>
void permute_rows(T *data, const CBLAS_INT *perm, T *tmp, const npy_intp m, const npy_intp ncols)
{
    std::memcpy(tmp, data, (std::size_t)m * ncols * sizeof(T));
    for (npy_intp i = 0; i < m; i++) {
        if (perm[i] != i) {
            std::memcpy(data + i * ncols, tmp + perm[i] * ncols, ncols * sizeof(T));
        }
    }
}


/**
 * @brief LU decomposition of one m-by-n matrix.
 *
 * ?getrf factors the column-major @p f_buf in place; L and U are then written to the
 * zero-initialised row-major @p l_out and @p u_out, and the row permutation to @p perm.
 *
 * @param[in,out] f_buf      Column-major m-by-n input. Overwritten by ?getrf and, when
 *                           @p permute_l is set, reused as scratch by permute_rows.
 * @param[out]    l_out      Row-major m-by-min(m, n) unit lower triangular L, zero on entry.
 * @param[out]    u_out      Row-major min(m, n)-by-n upper triangular U, zero on entry.
 * @param[out]    ipiv       Scratch of min(m, n) elements for the pivot sequence.
 * @param[out]    perm       Row permutation of length m with P @ A = L @ U.
 * @param[in]     m          Number of rows.
 * @param[in]     n          Number of columns.
 * @param[in]     permute_l  If true, @p l_out receives P @ L instead of L.
 * @return LAPACK info: 0 success; > 0 U(info, info) is exactly zero but the factorization
 *         is complete; < 0 argument error.
 */
template<typename T>
CBLAS_INT lu_decompose(T *f_buf, T *l_out, T *u_out, CBLAS_INT *ipiv, CBLAS_INT *perm, const CBLAS_INT m, const CBLAS_INT n, const bool permute_l)
{
    const CBLAS_INT mn = m < n ? m : n;
    CBLAS_INT info = 0;

    lapack::getrf(m, n, f_buf, m, ipiv, &info);
    if (info < 0) { return info; }

    // index arithmetic in npy_intp: a slice may hold more than CBLAS_INT_MAX elements
    // U: the upper triangle and the right block of the column-major factor
    for (npy_intp i = 0; i < mn; i++) {
        for (npy_intp j = i; j < n; j++) {
            u_out[i * n + j] = f_buf[j * m + i];
        }
    }
    // L: the strictly lower triangle plus a unit diagonal
    for (npy_intp j = 0; j < mn; j++) {
        for (npy_intp i = j + 1; i < m; i++) {
            l_out[i * mn + j] = f_buf[j * m + i];
        }
    }
    for (npy_intp i = 0; i < mn; i++) { l_out[i * mn + i] = T(1); }

    ipiv_to_perm(ipiv, perm, m, mn);
    if (permute_l) { permute_rows(l_out, perm, f_buf, m, mn); }   // f_buf is free as scratch

    return info;
}


/// The type a determinant of T is computed and returned in: double precision for f32 and c64.
template<typename T>
using det_type_t = std::conditional_t<std::is_same_v<T, f32>, f64,
                   std::conditional_t<std::is_same_v<T, c64>, c128, T>>;


/**
 * @brief Determinant of one n-by-n matrix from its LU factorization.
 *
 * The product of the diagonal is accumulated and returned in det_type_t<T>, so that single
 * precision input does not overflow or underflow in the product.
 *
 * @param[in,out] f_buf  Column-major n-by-n input; overwritten by ?getrf.
 * @param[out]    ipiv   Scratch of n elements for the pivot sequence.
 * @param[in]     n      Matrix dimension.
 * @param[out]    info   LAPACK info: 0 success; > 0 exactly singular; < 0 argument error.
 * @return The determinant; zero when @p info is nonzero.
 */
template<typename T>
det_type_t<T> det_from_lu(T *f_buf, CBLAS_INT *ipiv, const CBLAS_INT n, CBLAS_INT *info)
{
    using acc_type = det_type_t<T>;
    lapack::getrf(n, n, f_buf, n, ipiv, info);
    if (*info != 0) { return acc_type(0); }

    acc_type det = acc_type(1);
    CBLAS_INT swaps = 0;
    for (npy_intp k = 0; k < n; k++) {
        det *= static_cast<acc_type>(f_buf[k * (n + 1)]);  // diagonal entry (k, k), column-major
        if (ipiv[k] != k + 1) { swaps++; }
    }
    return (swaps % 2) ? -det : det;
}


/**
 * @brief Batched LU decomposition, `linalg.lu`.
 *
 * @param[in]  ap_a       Input array of shape (..., m, n), any strides.
 * @param[out] l_out      C-contiguous (..., m, min(m, n)) buffer for L, zero on entry.
 * @param[out] u_out      C-contiguous (..., min(m, n), n) buffer for U, zero on entry.
 * @param[out] perm       C-contiguous (..., m) buffer for the row permutations.
 * @param[in]  permute_l  If true, L receives P @ L; the permutation is still written.
 * @param[in]  in_place   If true, @p ap_a is 2-D, F-contiguous and writeable, and is factored
 *                        in place instead of being copied to a workspace.
 * @return ok; too_large if m or n does not fit CBLAS_INT; lapack_error if ?getrf rejects an
 *         argument. An exactly singular slice is not an error.
 * @throws std::bad_alloc if a workspace cannot be allocated.
 */
template<typename T>
BatchStatus
_lu(PyArrayObject *ap_a, T *l_out, T *u_out, CBLAS_INT *perm, const bool permute_l, const bool in_place)
{
    char *a_data = PyArray_BYTES(ap_a);
    npy_intp ndim = PyArray_NDIM(ap_a);
    npy_intp *shape = PyArray_SHAPE(ap_a);
    npy_intp *strides = PyArray_STRIDES(ap_a);
    npy_intp s0 = strides[ndim - 2], s1 = strides[ndim - 1];

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    CBLAS_INT m, n;
    if (!checked_narrow<CBLAS_INT>(shape[ndim - 2], &m) || !checked_narrow<CBLAS_INT>(shape[ndim - 1], &n)) {
        return BatchStatus::too_large;
    }
    const npy_intp mn = m < n ? m : n;

    // the count is bounded by an existing array, so it is exact in 64 bits
    auto scratch = alloc_workspace<T>(in_place ? 0 : (std::uint64_t)m * n);
    auto ipiv = alloc_workspace<CBLAS_INT>(mn);

    for (npy_intp idx = 0; idx < outer_size; idx++) {
        T *slice_ptr = compute_slice_ptr<T>(idx, a_data, ndim, shape, strides);
        T *f_buf = slice_ptr;
        if (!in_place) {
            copy_strided_to_f(scratch.get(), slice_ptr, m, n, s0, s1, m);
            f_buf = scratch.get();
        }

        CBLAS_INT info = lu_decompose(f_buf, l_out + idx * m * mn, u_out + idx * mn * n,
                                      ipiv.get(), perm + idx * m, m, n, permute_l);
        if (info < 0) { return BatchStatus::lapack_error; }
    }

    return BatchStatus::ok;
}


/**
 * @brief Batched determinant, `linalg.det`.
 *
 * @param[in]  ap_a      Input array of shape (..., n, n), any strides.
 * @param[out] det_out   C-contiguous (...) buffer, one determinant per slice, in det_type_t<T>.
 * @param[in]  in_place  If true, @p ap_a is 2-D, contiguous in either order and writeable,
 *                       and is factored in place; det(A^T) = det(A).
 * @return ok; too_large if n does not fit CBLAS_INT; lapack_error if ?getrf rejects an
 *         argument. An exactly singular slice has determinant zero.
 * @throws std::bad_alloc if a workspace cannot be allocated.
 */
template<typename T>
BatchStatus
_det(PyArrayObject *ap_a, det_type_t<T> *det_out, const bool in_place)
{
    char *a_data = PyArray_BYTES(ap_a);
    npy_intp ndim = PyArray_NDIM(ap_a);
    npy_intp *shape = PyArray_SHAPE(ap_a);
    npy_intp *strides = PyArray_STRIDES(ap_a);
    npy_intp s0 = strides[ndim - 2], s1 = strides[ndim - 1];

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    CBLAS_INT n;
    if (!checked_narrow<CBLAS_INT>(shape[ndim - 1], &n)) { return BatchStatus::too_large; }

    // the count is bounded by an existing array, so it is exact in 64 bits
    auto scratch = alloc_workspace<T>(in_place ? 0 : (std::uint64_t)n * n);
    auto ipiv = alloc_workspace<CBLAS_INT>(n);

    for (npy_intp idx = 0; idx < outer_size; idx++) {
        T *slice_ptr = compute_slice_ptr<T>(idx, a_data, ndim, shape, strides);
        T *f_buf = slice_ptr;
        if (!in_place) {
            copy_strided_to_f(scratch.get(), slice_ptr, n, n, s0, s1, n);
            f_buf = scratch.get();
        }

        CBLAS_INT info = 0;
        det_out[idx] = det_from_lu(f_buf, ipiv.get(), n, &info);
        if (info < 0) { return BatchStatus::lapack_error; }
    }

    return BatchStatus::ok;
}

} // namespace sp_linalg
