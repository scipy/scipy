/*
 * Templated loops for `linalg.svd`
 */
#pragma once

namespace sp_linalg {

/**
 * @brief Shapes of U and Vh for an m-by-n input: U is (m, m) or (m, k), Vh is (n, n) or (k, n).
 *
 * @param m          Number of rows of the input.
 * @param n          Number of columns of the input.
 * @param jobz       'N' (no U, Vh), 'A' (full U, Vh) or 'S' (economic U, Vh).
 * @param u_shape0   Receives the number of rows of U, 0 for 'N'.
 * @param u_shape1   Receives the number of columns of U, 0 for 'N'.
 * @param vh_shape0  Receives the number of rows of Vh, 0 for 'N'.
 * @param vh_shape1  Receives the number of columns of Vh, 0 for 'N'.
 */
inline void
u_vh_shapes(npy_intp m, npy_intp n, char jobz,
            npy_intp *u_shape0, npy_intp *u_shape1,
            npy_intp *vh_shape0, npy_intp *vh_shape1
){
    npy_intp k = m < n ? m : n;

    switch(jobz) {
        case('N') :
            *u_shape0 = 0;
            *u_shape1 = 0;
            *vh_shape0 = 0;
            *vh_shape1 = 0;
            break;
        case('A'):
            *u_shape0 = m;
            *u_shape1 = m;
            *vh_shape0 = n;
            *vh_shape1 = n;
            break;
        default:  // 'S'
            *u_shape0 = m;
            *u_shape1 = k;
            *vh_shape0 = k;
            *vh_shape1 = n;
            break;
    }
}


/**
 * @brief Minimum lwork of real ?gesdd for an m-by-n input, on the path ?gesdd takes for it.
 *
 * Mirrors MINWRK in LAPACK's sgesdd/dgesdd. A ?gesdd workspace query reports at least this value.
 *
 * @tparam real_type  f32 or f64; the precision ?gesdd evaluates its path threshold in.
 * @param m     Number of rows.
 * @param n     Number of columns.
 * @param jobz  'N', 'S' or 'A'.
 * @return The minimum lwork, evaluated without overflow.
 */
template<typename real_type>
std::uint64_t
gesdd_min_lwork(std::uint64_t m, std::uint64_t n, char jobz)
{
    const std::uint64_t k = std::min(m, n), kmax = std::max(m, n);
    const std::uint64_t bdspac = (jobz == 'N') ? 7*k : 3*k*k + 4*k;
    // the paths that start with a QR or LQ factorization (MNTHR in ?gesdd)
    const bool tall = kmax >= (std::uint64_t)(real_type(k) * real_type(11) / real_type(6));
    if (!tall) { return 3*k + std::max(kmax, bdspac); }
    if (jobz == 'N') { return bdspac + k; }
    if (jobz == 'S') { return bdspac + k*k + 3*k; }
    return k*k + std::max(3*k + bdspac, k + kmax);
}


template<typename T>
BatchStatus
_svd(PyArrayObject* ap_Am, PyArrayObject *ap_U, PyArrayObject *ap_S, PyArrayObject *ap_Vh, char jobz, bool use_gesdd, int overwrite_a, SliceStatusVec& vec_status)
{
    using real_type = real_of_t<T>; // f32 if T==c64 etc
    SliceStatus slice_status;

    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    char *Am_data = PyArray_BYTES(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp m = shape[ndim - 2];                // Slice size
    npy_intp n = shape[ndim - 1];
    npy_intp* strides = PyArray_STRIDES(ap_Am);

    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2) {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    // output array attributes
    npy_intp min_mn = m < n ? m : n;
    npy_intp max_mn = m < n ? n : m;

    // U, Vh shapes ( =0 if compute_uv=False)
    npy_intp u_shape0, u_shape1, vh_shape0, vh_shape1;
    u_vh_shapes(m, n, jobz, &u_shape0, &u_shape1, &vh_shape0, &vh_shape1);

    CBLAS_INT intn, intm, ldu, ldvh, info;
    if (!checked_narrow<CBLAS_INT>(n, &intn) || !checked_narrow<CBLAS_INT>(m, &intm) ||
        !checked_narrow<CBLAS_INT>(u_shape0 == 0 ? 1 : u_shape0, &ldu) ||
        !checked_narrow<CBLAS_INT>(vh_shape0 == 0 ? 1 : vh_shape0, &ldvh)) {
        return BatchStatus::too_large;
    }

    T *ptr_U = u_shape0 > 0 ? (T *)PyArray_DATA(ap_U) : NULL;
    T *ptr_Vh = vh_shape0 > 0 ? (T *)PyArray_DATA(ap_Vh) : NULL;
    real_type *ptr_S = (real_type *)PyArray_DATA(ap_S);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    // gesdd needs iwork of 8*min(m, n); the complex flavors need rwork, quadratic in min(m, n) for
    // gesdd with vectors
    const std::uint64_t k = min_mn, kmax = max_mn;
    std::uint64_t iwork_size = 0, rwork_size = 0;
    if (use_gesdd) {
        iwork_size = 8 * k;
        // ?gesdd evaluates its workspace in CBLAS_INT, already in the workspace query; with
        // vectors the largest of those sizes is the minimum of its path, quadratic in min(m, n)
        if constexpr (!is_complex_v<T>) {
            if (gesdd_min_lwork<real_type>(m, n, jobz) > (std::uint64_t)CBLAS_INT_MAX) { return BatchStatus::too_large; }
        }
        if constexpr (is_complex_v<T>) {
            // without vectors LAPACK <= 3.6 needs 7*min(m, n), later versions 5*min(m, n)
            rwork_size = (jobz == 'N') ? 7 * k : std::max(5*k*k + 5*k, 2*kmax*k + 2*k*k + k);
            // ?gesdd addresses rwork with CBLAS_INT offsets
            if (rwork_size > (std::uint64_t)CBLAS_INT_MAX) { return BatchStatus::too_large; }
        }
    } else {
        if constexpr (is_complex_v<T>) { rwork_size = 5 * k; }
    }

    T tmp = 0.0;

    // query LWORK
    if (use_gesdd) {
        if constexpr (!is_complex_v<T>) {
            lapack::gesdd(jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, NULL, &info);
        } else {
            lapack::gesdd(jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, NULL, NULL, &info);
        }
    } else {
        if constexpr (!is_complex_v<T>) {
            lapack::gesvd(jobz, jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, &info);
        } else {
            lapack::gesvd(jobz, jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, NULL, &info);
        }
    }
    if (info != 0) { return BatchStatus::lapack_error; }

    // _calc_lwork returns -1 when the queried lwork does not fit CBLAS_INT
    CBLAS_INT lwork = _calc_lwork(tmp);
    if (lwork < 0) { return BatchStatus::too_large; }

    /*
     * Allocate memory and chop the buffer into parts
     *
     *    lwork     data_size
     * |----------|-----------|----------|------|
     * ^          ^           ^          ^
     * work       data        buf_U     buf_Vh
     *
     * Here
     *   - data is m*n if overwrite_a is False (and =0 otherwise)
     *   - buf_U is u_shape0*u_shape1 if jobz !='N' (and =0 otherwise)
     *   - buf_Vh is vh_shape0*vh_shapw1 if jobz != 'N' (and =0 otherwise)
     *
     * NB: we do not implement jobz='O' yet, so we never reuse A for U or Vh.
     */
    // each count is bounded by an existing array (A, U, Vh) or by CBLAS_INT, so the 64-bit sums are exact
    const std::uint64_t data_size = overwrite_a ? 0 : (std::uint64_t)m * n;
    const std::uint64_t u_size = (std::uint64_t)u_shape0 * u_shape1, vh_size = (std::uint64_t)vh_shape0 * vh_shape1;

    auto buf = alloc_workspace<T>(lwork + data_size + u_size + vh_size);
    auto iwork = alloc_workspace<CBLAS_INT>(iwork_size);
    auto rwork = alloc_workspace<real_type>(rwork_size);

    // partition the workspace
    T *cursor = buf.get();
    T *work = cursor; cursor += lwork;
    T *data = static_cast<T *>(PyArray_DATA(ap_Am));
    if (!overwrite_a) { data = cursor; cursor += data_size; }
    T *buf_U = NULL, *buf_Vh = NULL;
    if (jobz != 'N') { buf_U = cursor; cursor += u_size; buf_Vh = cursor; }

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            // copy the slice to `data` in F order
            T *slice_ptr = compute_slice_ptr<T>(idx, Am_data, ndim, shape, strides);
            copy_strided_to_f(data, slice_ptr, m, n, strides[ndim-2], strides[ndim-1], m);
        }

        // SVD the slice
        if (use_gesdd) {
            if constexpr (!is_complex_v<T>) {
                lapack::gesdd(jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, iwork.get(), &info);
            } else {
                lapack::gesdd(jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, rwork.get(), iwork.get(), &info);
            }
        } else {
            if constexpr (!is_complex_v<T>) {
                lapack::gesvd(jobz, jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, &info);
            } else {
                lapack::gesvd(jobz, jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, rwork.get(), &info);
            }
        }

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            return BatchStatus::ok;
        }

        // copy-and-transpose U and Vh slices from temp buffers to the output;
        // Also advance the output pointers: U, S, Vh are C-contiguous by construction
        if (jobz != 'N') {
            copy_slice_F_to_C(ptr_U, buf_U, u_shape0, u_shape1);
            ptr_U += u_shape0 * u_shape1;

            copy_slice_F_to_C(ptr_Vh, buf_Vh, vh_shape0, vh_shape1);
            ptr_Vh += vh_shape0 * vh_shape1;
        }

        // S slice has been filled in-place already
        ptr_S += min_mn;
    }

    return BatchStatus::ok;
}

} // namespace sp_linalg
