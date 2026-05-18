#pragma once
#include <cstring>
#include <cstdint>
#include <type_traits>
#include "scipy_blas_defines.h"

namespace sp_linalg {


constexpr int LU_MAX_NDIM = 64;


// ========================================================================
// LAPACK dispatch
// ========================================================================

/**
 * @brief Dispatch to the appropriate LAPACK ?getrf routine based on type.
 *
 * This header does not declare its own extern "C" getrf symbols — it reuses
 * the declarations from _common_array_utils.hh (which must be included first).
 * Those declarations use npy_complex types, so std::complex pointers are
 * reinterpret_cast'd here; the types are ABI-compatible.
 *
 * We provide our own template wrapper (rather than using call_getrf from
 * _common_array_utils.hh) because call_getrf's overloads take npy_complex
 * parameters, and this header works exclusively with std::complex types
 * to keep NumPy dependencies out of the C++ logic.
 *
 * @param[in]     m     Number of rows.
 * @param[in]     n     Number of columns.
 * @param[in,out] a     Column-major m-by-n matrix; overwritten with L and U on exit.
 * @param[in]     lda   Leading dimension of @p a (>= m).
 * @param[out]    ipiv  Pivot indices, length min(m, n). Row i was interchanged with row ipiv[i].
 * @param[out]    info  0 on success, < 0 for illegal argument, > 0 if U(info,info) is zero.
 */
template<typename T>
inline void getrf(CBLAS_INT *m, CBLAS_INT *n, T *a, CBLAS_INT *lda, CBLAS_INT *ipiv, CBLAS_INT *info)
{
    if      constexpr (std::is_same_v<T, float>)                BLAS_FUNC(sgetrf)(m, n, a, lda, ipiv, info);
    else if constexpr (std::is_same_v<T, double>)               BLAS_FUNC(dgetrf)(m, n, a, lda, ipiv, info);
    else if constexpr (std::is_same_v<T, std::complex<float>>)  BLAS_FUNC(cgetrf)(m, n, reinterpret_cast<npy_complex64*>(a), lda, ipiv, info);
    else if constexpr (std::is_same_v<T, std::complex<double>>) BLAS_FUNC(zgetrf)(m, n, reinterpret_cast<npy_complex128*>(a), lda, ipiv, info);
}


// ========================================================================
// Context struct
// ========================================================================

/**
 * @brief Metadata for batched LU decomposition.
 *
 * Filled by the Python C extension from PyArrayObject metadata.
 * No NumPy types or Python objects cross this boundary; all data
 * pointers are passed separately as template-typed arguments.
 *
 * @note Element strides (not byte strides) are stored in @c strides.
 */
struct LU_Context {
    int64_t*    shape;            ///< Full nD shape of the input array.
    int64_t*    strides;          ///< Element strides for the input array.
    int64_t     num_of_slices;    ///< Product of batch dimensions (shape[0..ndim-3]).
    CBLAS_INT   m, n;             ///< Rows and columns of each 2D slice.
    CBLAS_INT   *perm;            ///< Output permutation buffer, num_of_slices * m elements.
    CBLAS_INT   *ipiv;            ///< Pivot scratch buffer, min(m,n) elements (reused across slices).
    int         ndim;             ///< Number of dimensions in the input array.
    bool        permute_l;        ///< If true, apply P to L and discard the permutation.
    bool        overwrite_a;      ///< If true, input buffer is used directly as getrf workspace.
};


// ========================================================================
// Utility functions
// ========================================================================

/**
 * @brief Copy a strided m-by-n matrix into a contiguous column-major buffer.
 *
 * @param[in]  src         Source data pointer.
 * @param[out] dst         Destination buffer, at least m*n elements.
 * @param[in]  m           Number of rows.
 * @param[in]  n           Number of columns.
 * @param[in]  stride_row  Row stride in elements (not bytes).
 * @param[in]  stride_col  Column stride in elements (not bytes).
 */
template<typename T>
void copy_strided_to_f(const T *src, T *dst, CBLAS_INT m, CBLAS_INT n, CBLAS_INT stride_row, CBLAS_INT stride_col)
{
    for (CBLAS_INT i = 0; i < m; i++)
        for (CBLAS_INT j = 0; j < n; j++)
            dst[j * m + i] = src[i * stride_row + j * stride_col];
}


/**
 * @brief Convert LAPACK 1-based pivot swap sequence to a 0-based permutation.
 *
 * On exit, @p perm contains the permutation such that P @ A = L @ U.
 *
 * @param[in,out] ipiv  On entry, length-mn pivot sequence from ?getrf (1-based).
 *                      Used as scratch of length m during inversion. Contents are destroyed.
 * @param[out]    perm  Output permutation array, length m.
 * @param[in]     m     Number of rows.
 * @param[in]     mn    min(m, n) — length of the pivot sequence.
 */
inline void ipiv_to_perm(CBLAS_INT *ipiv, CBLAS_INT *perm, CBLAS_INT m, CBLAS_INT mn)
{
    // Apply the swap sequence in reverse to get the inverse permutation
    // directly, without extra scratch. Forward swaps give P, reverse gives P⁻¹.
    for (CBLAS_INT i = 0; i < m; i++) perm[i] = i;
    for (CBLAS_INT i = mn - 1; i >= 0; i--) {
        CBLAS_INT tmp = perm[i];
        perm[i] = perm[ipiv[i] - 1];
        perm[ipiv[i] - 1] = tmp;
    }
}


/**
 * @brief Permute rows of a C-ordered m-by-ncols matrix in place.
 *
 * On exit, row i of @p data contains original row perm[i].
 *
 * @param[in,out] data   Row-major matrix, m * ncols elements.
 * @param[in]     perm   Permutation array, length m.
 * @param[in]     tmp    Scratch buffer, at least m * ncols elements.
 * @param[in]     m      Number of rows.
 * @param[in]     ncols  Number of columns.
 */
template<typename T>
void permute_rows(T *data, const CBLAS_INT *perm, T *tmp, CBLAS_INT m, CBLAS_INT ncols)
{
    std::memcpy(tmp, data, m * ncols * sizeof(T));
    for (CBLAS_INT i = 0; i < m; i++) {
        if (perm[i] != i) {
            std::memcpy(data + i * ncols, tmp + perm[i] * ncols, ncols * sizeof(T));
        }
    }
}


// ========================================================================
// Single-slice LU decomposition
// ========================================================================

/**
 * @brief LU decomposition of a single m-by-n matrix.
 *
 * Calls ?getrf on the column-major buffer @p f_buf, then extracts L and U
 * into pre-zeroed row-major output buffers. Optionally permutes L rows so
 * that the caller receives P @ L directly.
 *
 * All buffers are caller-allocated and must satisfy the documented sizes.
 *
 * @param[in,out] f_buf      Column-major m*n buffer. On entry, the input matrix.
 *                           Overwritten by ?getrf. Reused as scratch for permute_rows.
 * @param[out]    l_out      Row-major m*min(m,n) buffer, pre-zeroed by caller.
 * @param[out]    u_out      Row-major min(m,n)*n buffer, pre-zeroed by caller.
 * @param[in,out] ipiv       Scratch buffer, at least min(m,n) elements. Contents destroyed.
 * @param[out]    perm       Output permutation, m elements.
 * @param[in]     m          Number of rows.
 * @param[in]     n          Number of columns.
 * @param[in]     permute_l  If true, apply P to L in place.
 *
 * @return LAPACK info: 0 = success, < 0 = bad argument, > 0 = U(info,info) is zero (singular).
 */
template<typename T>
CBLAS_INT lu_decompose(T *f_buf, T *l_out, T *u_out, CBLAS_INT *ipiv, CBLAS_INT *perm, CBLAS_INT m, CBLAS_INT n, bool permute_l)
{
    CBLAS_INT mn = m < n ? m : n;
    CBLAS_INT info = 0;

    getrf(&m, &n, f_buf, &m, ipiv, &info);

    if (info < 0) { return info; }

    // Extract U (upper triangle + right block) from F-ordered getrf result
    // into pre-zeroed C-ordered output buffer.
    for (CBLAS_INT i = 0; i < mn; i++) {
        for (CBLAS_INT j = i; j < n; j++) {
            u_out[i * n + j] = f_buf[j * m + i];
        }
    }

    // Extract L (below-diagonal elements) and set unit diagonal.
    for (CBLAS_INT j = 0; j < mn; j++)
        for (CBLAS_INT i = j + 1; i < m; i++)
            l_out[i * mn + j] = f_buf[j * m + i];

    for (CBLAS_INT i = 0; i < mn; i++)
        l_out[i * mn + i] = T(1);

    // Convert LAPACK pivot swaps to permutation
    ipiv_to_perm(ipiv, perm, m, mn);

    // Optionally permute L rows (f_buf is free to use as scratch)
    if (permute_l)
        permute_rows(l_out, perm, f_buf, m, mn);

    return info;
}


// ========================================================================
// Batched dispatch
// ========================================================================

/**
 * @brief Batched LU decomposition over an nD array.
 *
 * Iterates over all 2D slices defined by the batch dimensions,
 * copies each strided input slice into the column-major scratch buffer,
 * and calls lu_decompose(). Per-slice LAPACK info is written to
 * @p slice_info for the Python layer to convert to warnings.
 *
 * When @c ctx.overwrite_a is true (2D, F-contiguous input only), the
 * input data is used directly as the ?getrf workspace — no copy is made.
 *
 * @param[in]     ctx         Batch metadata (shape, strides, dimensions, flags).
 * @param[in]     a_dat       Input array data, possibly strided.
 * @param[out]    l_out       Contiguous output buffer for L factors, num_of_slices * m * min(m,n).
 * @param[out]    u_out       Contiguous output buffer for U factors, num_of_slices * min(m,n) * n.
 * @param[in,out] scratch     Column-major work buffer, m * n elements. Unused when overwrite_a is true.
 * @param[out]    slice_info  Per-slice LAPACK info, num_of_slices elements. 0 = success, > 0 = singular.
 *
 * @return 0 on success, < 0 on fatal LAPACK error (processing stops immediately).
 */
template<typename T>
int lu_dispatch(LU_Context &ctx, T *a_dat, T *l_out, T *u_out, T *scratch, CBLAS_INT *slice_info)
{
    CBLAS_INT m = ctx.m, n = ctx.n;
    CBLAS_INT mn = m < n ? m : n;

    CBLAS_INT s0 = (CBLAS_INT)ctx.strides[ctx.ndim - 2];
    CBLAS_INT s1 = (CBLAS_INT)ctx.strides[ctx.ndim - 1];

    for (int64_t idx = 0; idx < ctx.num_of_slices; idx++) {

        // Compute strided offset for this slice's input data
        int64_t offset = 0;
        int64_t temp = idx;
        for (int i = ctx.ndim - 3; i >= 0; i--) {
            offset += (temp % ctx.shape[i]) * ctx.strides[i];
            temp /= ctx.shape[i];
        }
        const T *slice_in = a_dat + offset;

        // Contiguous output offsets
        T *l_slice = l_out + idx * m * mn;
        T *u_slice = u_out + idx * mn * n;
        CBLAS_INT *perm_slice = ctx.perm + idx * m;

        // Prepare column-major input for getrf
        T *f_buf;
        if (ctx.overwrite_a) {
            // 2D F-contiguous: use input directly (no copy)
            f_buf = const_cast<T*>(slice_in);
        } else {
            copy_strided_to_f(slice_in, scratch, m, n, s0, s1);
            f_buf = scratch;
        }

        // Decompose this slice
        CBLAS_INT info = lu_decompose(f_buf, l_slice, u_slice,
                                      ctx.ipiv, perm_slice, m, n,
                                      ctx.permute_l);

        if (info < 0) { return (int)info; }
        slice_info[idx] = info;
    }

    return 0;
}


// ========================================================================
// Single-slice determinant
// ========================================================================

/**
 * @brief Compute the determinant of a single n-by-n matrix via LU factorization.
 *
 * Calls ?getrf on the column-major buffer @p f_buf, then computes the
 * product of the diagonal of U and adjusts the sign based on the
 * number of row interchanges.
 *
 * @param[in,out] f_buf  Column-major n*n buffer. Overwritten by ?getrf.
 * @param[in,out] ipiv   Scratch buffer, at least n elements. Contents destroyed.
 * @param[in]     n      Matrix dimension.
 * @param[out]    info   LAPACK info: 0 = success, < 0 = bad arg, > 0 = singular.
 *
 * @return The determinant value. Returns T(0) when singular (info > 0).
 */
template<typename T>
T det_from_lu(T *f_buf, CBLAS_INT *ipiv, CBLAS_INT n, CBLAS_INT *info)
{
    getrf(&n, &n, f_buf, &n, ipiv, info);

    if (*info < 0) { return T(0); }
    if (*info > 0) { return T(0); }

    // Accumulate in promoted precision to avoid overflow/underflow
    // for single-precision types (float32 -> float64, complex64 -> complex128).
    using acc_type = std::conditional_t<
        std::is_same_v<T, float>,
        double,
        std::conditional_t<
            std::is_same_v<T,std::complex<float>>,
            std::complex<double>,
            T
        >
    >;

    acc_type det = acc_type(1);
    CBLAS_INT swaps = 0;
    for (CBLAS_INT k = 0; k < n; k++) {
        det *= static_cast<acc_type>(f_buf[k * (n + 1)]);  // diagonal (k,k) in column-major
        if (ipiv[k] != k + 1) { swaps++; }
    }
    return static_cast<T>((swaps % 2) ? -det : det);
}


// ========================================================================
// Batched determinant dispatch
// ========================================================================

/**
 * @brief Batched determinant computation over an nD array.
 *
 * Iterates over all 2D slices defined by the batch dimensions,
 * copies each strided input slice into the column-major scratch buffer,
 * and calls det_from_lu().
 *
 * When @c ctx.overwrite_a is true (2D, F-contiguous input only), the
 * input data is used directly as the ?getrf workspace.
 *
 * @param[in]     ctx         Batch metadata (shape, strides, dimensions, flags).
 * @param[in]     a_dat       Input array data, possibly strided.
 * @param[out]    det_out     Contiguous output buffer, one element per slice.
 * @param[in,out] scratch     Column-major work buffer, n * n elements. Unused when overwrite_a.
 * @param[out]    slice_info  Per-slice LAPACK info, num_of_slices elements.
 *
 * @return 0 on success, < 0 on fatal LAPACK error (processing stops immediately).
 */
template<typename T>
int det_dispatch(LU_Context &ctx, T *a_dat, T *det_out, T *scratch, CBLAS_INT *slice_info)
{
    CBLAS_INT n = ctx.n;

    CBLAS_INT s0 = (CBLAS_INT)ctx.strides[ctx.ndim - 2];
    CBLAS_INT s1 = (CBLAS_INT)ctx.strides[ctx.ndim - 1];

    for (int64_t idx = 0; idx < ctx.num_of_slices; idx++) {

        // Compute strided offset for this slice's input data
        int64_t offset = 0;
        int64_t temp = idx;
        for (int i = ctx.ndim - 3; i >= 0; i--) {
            offset += (temp % ctx.shape[i]) * ctx.strides[i];
            temp /= ctx.shape[i];
        }
        const T *slice_in = a_dat + offset;

        // Prepare column-major input for getrf
        T *f_buf;
        if (ctx.overwrite_a) {
            f_buf = const_cast<T*>(slice_in);
        } else {
            copy_strided_to_f(slice_in, scratch, n, n, s0, s1);
            f_buf = scratch;
        }

        CBLAS_INT info = 0;
        det_out[idx] = det_from_lu(f_buf, ctx.ipiv, n, &info);

        if (info < 0) { return (int)info; }
        slice_info[idx] = info;
    }

    return 0;
}


} // namespace sp_linalg
