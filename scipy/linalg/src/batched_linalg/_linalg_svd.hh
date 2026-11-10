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
_svd_gesdd(PyArrayObject* ap_Am, PyArrayObject *ap_U, PyArrayObject *ap_S, PyArrayObject *ap_Vh, char jobz, int overwrite_a, SliceStatusVec& vec_status)
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
    if constexpr (!is_complex_v<T>) {
        lapack::gesdd(jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, NULL, &info);
    } else {
        lapack::gesdd(jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, NULL, NULL, &info);
    }
    if (info != 0) { return BatchStatus::lapack_query_failed; }

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
    std::size_t data_size = 0, u_size, vh_size, bufsize, iwork_size, rwork_size = 0;
    if ((!overwrite_a && !checked_mul(m, n, &data_size)) ||
        !checked_mul(u_shape0, u_shape1, &u_size) || !checked_mul(vh_shape0, vh_shape1, &vh_size) ||
        !checked_add(data_size, lwork, &bufsize) || !checked_add(bufsize, u_size, &bufsize) ||
        !checked_add(bufsize, vh_size, &bufsize) || !checked_mul(8, min_mn, &iwork_size)) {
        return BatchStatus::too_large;
    }
    if constexpr (is_complex_v<T>) {
        // assume LAPACK > 3.6 (cf LAPACK docs on netlib.org):
        // lrwork = max(5*min_mn*min_mn + 5*min_mn, 2*max_mn*min_mn + 2*min_mn*min_mn + min_mn)
        std::size_t kk, mk, five_k, two_kk, lrwork1, lrwork2;
        if (!checked_mul(min_mn, min_mn, &kk) || !checked_mul(max_mn, min_mn, &mk) ||
            !checked_mul(5, min_mn, &five_k) || !checked_mul(2, kk, &two_kk) ||
            !checked_mul(5, kk, &lrwork1) || !checked_add(lrwork1, five_k, &lrwork1) ||
            !checked_mul(2, mk, &lrwork2) || !checked_add(lrwork2, two_kk, &lrwork2) ||
            !checked_add(lrwork2, min_mn, &lrwork2)) {
            return BatchStatus::too_large;
        }
        rwork_size = std::max(lrwork1, lrwork2);
    }

    auto buf = alloc_workspace<T>(bufsize);
    auto iwork = alloc_workspace<CBLAS_INT>(iwork_size);
    auto rwork = alloc_workspace<real_type>(rwork_size);
    if (!buf || !iwork || !rwork) { return BatchStatus::no_memory; }

    // partition the workspace
    T *work = buf.get();
    T *data = overwrite_a ? Am_data : work + lwork;
    T *buf_U = (jobz != 'N') ? work + lwork + data_size : NULL;
    T *buf_Vh = (jobz != 'N') ? buf_U + u_size : NULL;

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            // copy the slice to `data` in F order
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_strided_to_f(slice_ptr, data, m, n, strides[ndim-2], strides[ndim-1], m);
        }

        // SVD the slice
        if constexpr (!is_complex_v<T>) {
            lapack::gesdd(jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, iwork.get(), &info);
        } else {
            lapack::gesdd(jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, rwork.get(), iwork.get(), &info);
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


template<typename T>
BatchStatus
_svd_gesvd(PyArrayObject* ap_Am, PyArrayObject *ap_U, PyArrayObject *ap_S, PyArrayObject *ap_Vh, char jobz, int overwrite_a, SliceStatusVec& vec_status)
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
    if constexpr (!is_complex_v<T>) {
        lapack::gesvd(jobz, jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, &info);
    } else {
        lapack::gesvd(jobz, jobz, intm, intn, NULL, intm, NULL, NULL, ldu, NULL, ldvh, &tmp, -1, NULL, &info);
    }
    if (info != 0) { return BatchStatus::lapack_query_failed; }

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
    std::size_t data_size = 0, u_size, vh_size, bufsize, rwork_size = 0;
    if ((!overwrite_a && !checked_mul(m, n, &data_size)) ||
        !checked_mul(u_shape0, u_shape1, &u_size) || !checked_mul(vh_shape0, vh_shape1, &vh_size) ||
        !checked_add(data_size, lwork, &bufsize) || !checked_add(bufsize, u_size, &bufsize) ||
        !checked_add(bufsize, vh_size, &bufsize) ||
        (is_complex_v<T> && !checked_mul(5, min_mn, &rwork_size))) {
        return BatchStatus::too_large;
    }

    auto buf = alloc_workspace<T>(bufsize);
    auto rwork = alloc_workspace<real_type>(rwork_size);
    if (!buf || !rwork) { return BatchStatus::no_memory; }

    // partition the workspace
    T *work = buf.get();
    T *data = overwrite_a ? Am_data : work + lwork;
    T *buf_U = (jobz != 'N') ? work + lwork + data_size : NULL;
    T *buf_Vh = (jobz != 'N') ? buf_U + u_size : NULL;

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            // copy the slice to `data` in F order
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_strided_to_f(slice_ptr, data, m, n, strides[ndim-2], strides[ndim-1], m);
        }

        // SVD the slice
        if constexpr (!is_complex_v<T>) {
            lapack::gesvd(jobz, jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, &info);
        } else {
            lapack::gesvd(jobz, jobz, intm, intn, data, intm, ptr_S, buf_U, ldu, buf_Vh, ldvh, work, lwork, rwork.get(), &info);
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

        // S slice has been filled in-place already;
        ptr_S += min_mn;
    }

    return BatchStatus::ok;
}

} // namespace sp_linalg
