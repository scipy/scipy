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


template<typename T>
BatchStatus
_svd(PyArrayObject* ap_Am, PyArrayObject *ap_U, PyArrayObject *ap_S, PyArrayObject *ap_Vh, char jobz, bool use_gesdd, int overwrite_a, SliceStatusVec& vec_status)
{
    using real_type = real_of_t<T>; // f32 if T==c64 etc
    SliceStatus slice_status;

    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    T* Am_data = (T *)PyArray_DATA(ap_Am);
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
    // element counts are products of dimensions of existing arrays (A, U, Vh), so they fit size_t
    const std::size_t data_size = overwrite_a ? 0 : (std::size_t)m * n;
    const std::size_t u_size = (std::size_t)u_shape0 * u_shape1, vh_size = (std::size_t)vh_shape0 * vh_shape1;

    // gesdd needs iwork of 8*min_mn; complex flavors need rwork, much larger for gesdd
    std::size_t iwork_size = 0, rwork_size = 0;
    if (use_gesdd) {
        iwork_size = 8 * (std::size_t)min_mn;
        if constexpr (is_complex_v<T>) {
            // assume LAPACK > 3.6 (cf LAPACK docs on netlib.org); at most about 5*m*n
            const std::size_t k = min_mn, kmax = max_mn;
            rwork_size = std::max(5*k*k + 5*k, 2*kmax*k + 2*k*k + k);
        }
    } else {
        if constexpr (is_complex_v<T>) { rwork_size = 5 * (std::size_t)min_mn; }
    }

    auto buf = alloc_workspace<T>(lwork + data_size + u_size + vh_size);
    auto iwork = alloc_workspace<CBLAS_INT>(iwork_size);
    auto rwork = alloc_workspace<real_type>(rwork_size);
    if (!buf || !iwork || !rwork) { return BatchStatus::no_memory; }

    // partition the workspace
    T *cursor = buf.get();
    T *work = cursor; cursor += lwork;
    T *data = Am_data;
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
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
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
