/*
 * Templated loops for `linalg.qr`
 */
#pragma once

namespace sp_linalg {

template<typename T>
BatchStatus
_qr(PyArrayObject *ap_Am, PyArrayObject *ap_Q, PyArrayObject *ap_R, PyArrayObject *ap_tau, PyArrayObject *ap_jpvt, bool overwrite_a, QR_mode mode, bool pivoting, SliceStatusVec &vec_status)
{
    using real_type = real_of_t<T>;
    SliceStatus slice_status;

    // ------------------------------------------------------------------------
    // Input array attributes, some bookkeeping due to conditional allocations
    // ------------------------------------------------------------------------
    char *A_data = PyArray_BYTES(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp M = shape[ndim-2];
    npy_intp N = shape[ndim-1];

    char *R_data = PyArray_BYTES(ap_R);
    npy_intp *strides_R = PyArray_STRIDES(ap_R);

    // Conditional allocations based on the mode; not all data is always available.
    char *Q_data = NULL, *tau_data = NULL, *jpvt_data = NULL;
    npy_intp *strides_Q = NULL, *strides_tau = NULL, *strides_jpvt = NULL;

    if (mode != QR_mode::R) {
        Q_data = PyArray_BYTES(ap_Q);
        strides_Q = PyArray_STRIDES(ap_Q);
    }

    if (mode == QR_mode::RAW_MODE) {
        tau_data = PyArray_BYTES(ap_tau);
        strides_tau = PyArray_STRIDES(ap_tau);
    }

    if (pivoting) {
        jpvt_data = PyArray_BYTES(ap_jpvt);
        strides_jpvt = PyArray_STRIDES(ap_jpvt);
    }

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    // -------------------------------------------------------------------
    // Workspace computation and allocation
    // -------------------------------------------------------------------
    CBLAS_INT intn, intm, info = 0;
    if (!checked_narrow<CBLAS_INT>(N, &intn) || !checked_narrow<CBLAS_INT>(M, &intm)) {
        return BatchStatus::too_large;
    }
    CBLAS_INT K = std::min(intn, intm), max_dim = std::max(intn, intm);
    CBLAS_INT middle_dim = (mode == QR_mode::ECONOMIC || mode == QR_mode::RAW_MODE) ? K : intm; // Final dimension: Q is [`M` x `middle_dim`], R is [`middle_dim` x `N`]


    // Probe both the factorization as well as `or_un_gqr` to find the optimal lwork
    T tmp_factor = 0.0;
    T tmp_or_un_gqr = 0.0;
    CBLAS_INT lwork = -1;

    if (!pivoting) {
        lapack::geqrf(intm, intn, NULL, intm, NULL, &tmp_factor, lwork, &info);
    } else {
        if constexpr (!is_complex_v<T>) {
            lapack::geqp3(intm, intn, NULL, intm, NULL, NULL, &tmp_factor, lwork, &info);
        } else {
            lapack::geqp3(intm, intn, NULL, intm, NULL, NULL, &tmp_factor, lwork, NULL, &info);
        }
    }
    if (info != 0) { return BatchStatus::lapack_error; }

    // Also probe the Q construction step if required.
    if (mode == QR_mode::FULL || mode == QR_mode::ECONOMIC) {
        if constexpr (!is_complex_v<T>) {
            lapack::orgqr(intm, middle_dim, K, NULL, intm, NULL, &tmp_or_un_gqr, lwork, &info);
        } else {
            lapack::ungqr(intm, middle_dim, K, NULL, intm, NULL, &tmp_or_un_gqr, lwork, &info);
        }
    }
    if (info != 0) { return BatchStatus::lapack_error; }

    // _calc_lwork returns -1 when the queried lwork does not fit CBLAS_INT
    CBLAS_INT lwork_factor = _calc_lwork(tmp_factor);
    CBLAS_INT lwork_or_un_gqr = _calc_lwork(tmp_or_un_gqr);
    if ((lwork_factor < 0) || (lwork_or_un_gqr < 0)) { return BatchStatus::too_large; }
    lwork = std::max(lwork_factor, lwork_or_un_gqr);

    // `geqp3` requires an `lwork` of at least 3N + 1, which must fit CBLAS_INT as well
    if (pivoting && (3 * (std::uint64_t)N + 1 > (std::uint64_t)CBLAS_INT_MAX)) { return BatchStatus::too_large; }

    /*
     * Allocate buffer and slice into parts
     *
     *     lwork        tau_size     buf_size
     * |------------|-------------|--------------|
     * ^            ^             ^
     * work         tau_buffer    data_A
     *
     * - `work` is the work area, size `lwork`. Always allocated.
     *
     * - `tau_buffer` contains a temporary buffer for the reflectors. If mode == RAW_MODE, they have
     *   to be returned and hence a temporary buffer is not needed. Else it has size `K`.
     *
     * - `data_A` contains a temporary buffer for LAPACK calls, if mode == FULL, the resulting
     *   `Q` is `M x max_dim` to account for both `M >= N` and `M < N`. For mode == R and
     *   mode == ECONOMIC the required buffer of size `M x N`. For mode == RAW_MODE this buffer
     *   can be bypassed altogether due to the fact that the returned output has F-ordered slices.
     *
     *   When `overwrite_a` is set the buffer is also not needed for ECONOMIC and R modes. Due to the fact
     *   that the resulting `Q` can be larger than `A` for mode == FULL a buffer is still allocated in that
     *   case.
     *
     * The re-use of the input buffer when `overwrite_a` is enabled works as follows:
     * - mode == RAW_MODE:
     *      the factorization is performed in-place and `R` is extracted into a newly allocated array.
     * - mode == R:
     *      factorization performed in-place and the other triangle is zeroed out.
     * - mode == ECONOMIC:
     *      perform factorization in-place, extract `R` to newly allocated array and explicitly compute
     *      `Q` in-place as well. If needed (i.e. `M < N`), the `Q` is shrunk into having `M x M` slices
     *      at the python side.
     * - mode == FULL:
     *      copy data into the temporary buffer and perform factorizations there. The input array is re-used for
     *      storing `R`. At the end `Q` is copied into a newly allocated array.
     *      XXX: for `M < N` it is possible to circumvent the temporary buffer, but that would lead to shape-dependent codepaths.
     */
    // each count is bounded by an existing array (A, and Q for mode == FULL) or by CBLAS_INT, so
    // the 64-bit sums are exact
    const std::uint64_t tau_size = (mode == QR_mode::RAW_MODE) ? 0 : K;
    const bool needs_buffer = !((overwrite_a && mode != QR_mode::FULL) || mode == QR_mode::RAW_MODE);
    const std::uint64_t buf_size = needs_buffer ? (std::uint64_t)M * ((mode == QR_mode::FULL) ? max_dim : N) : 0;
    const std::uint64_t rwork_size = (pivoting && is_complex_v<T>) ? 2 * (std::uint64_t)N : 0;

    auto buffer = alloc_workspace<T>(lwork + tau_size + buf_size);
    auto rwork = alloc_workspace<real_type>(rwork_size);   // `c/zgeqp3` needs rwork

    T *cursor = buffer.get();
    T *work = cursor;       cursor += lwork;
    T *tau_buffer = cursor; cursor += tau_size;
    T *data_A = cursor;

    // -------------------------------------------------------------------
    // Main loop to traverse the slices
    // A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
    // -------------------------------------------------------------------
    T *slice_ptr_A = NULL, *slice_ptr_Q = NULL, *slice_ptr_R = NULL, *slice_ptr_tau = NULL;
    CBLAS_INT *slice_ptr_jpvt = NULL;
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        // Compute slice pointers if necessary. `tau` and `jpvt` are vectors, so actually 1 dimension
        // less, however, the looping to compute the slices should traverse the same number of dimensions.
        // Hence, the use of `ndim` instead of `ndim-1`. Similarly, the shapes of `Q` and `R` might differ,
        // but only the batching ones are relevant, which are identical.
        slice_ptr_A = compute_slice_ptr<T>(idx, A_data, ndim, shape, strides);
        if (!(mode == QR_mode::R && overwrite_a)) {
            slice_ptr_R = compute_slice_ptr<T>(idx, R_data, ndim, shape, strides_R);
        } // NB. if the condition is true the relevant buffer is already in place and no copy back is needed either.

        if (mode != QR_mode::R) {
            slice_ptr_Q = compute_slice_ptr<T>(idx, Q_data, ndim, shape, strides_Q);
        }
        if (mode == QR_mode::RAW_MODE) {
            slice_ptr_tau = compute_slice_ptr<T>(idx, tau_data, ndim, shape, strides_tau);
        } else { // `tau` might still be required for computation of `Q`, but should not be returned, so the temporary buffer is sufficient.
            slice_ptr_tau = tau_buffer;
        }
        if (pivoting) {
            slice_ptr_jpvt = compute_slice_ptr<CBLAS_INT>(idx, jpvt_data, ndim, shape, strides_jpvt);
        }

        /*
         * Copy the buffer for processing. When `overwrite_a` is disabled the data should be handled in
         * another buffer, so copying is necessary. When the flag is set, the data can be processed in place.
         *
         * Another optimization is to process the data for mode == RAW in the return object directly as
         * the returned array is in F-order, hence bypassing the buffer avoids two redundant copies.
         */
        if ((!overwrite_a && mode != QR_mode::RAW_MODE) || mode == QR_mode::FULL) {
            copy_strided_to_f(data_A, slice_ptr_A, M, N, strides[ndim-2], strides[ndim-1], M);
        } else if (!overwrite_a && mode == QR_mode::RAW_MODE) {
            copy_strided_to_f(slice_ptr_Q, slice_ptr_A, M, N, strides[ndim-2], strides[ndim-1], M);
            data_A = slice_ptr_Q; // Add alias for easier to read codepaths later on.
        } else if (overwrite_a && (mode == QR_mode::RAW_MODE || mode == QR_mode::R || mode == QR_mode::ECONOMIC)) {
            data_A = slice_ptr_A;
        }

        init_status(slice_status, idx, St::GENERAL);

        // ---------------------------------------------------------------
        // Heavy lifting: call appropriate LAPACK routines.
        // ---------------------------------------------------------------

        // Factorization step, identical for each of the modes so do jointly.
        if (!pivoting) {
            lapack::geqrf(intm, intn, data_A, intm, slice_ptr_tau, work, lwork, &info);
        } else {
            if constexpr (!is_complex_v<T>) {
                lapack::geqp3(intm, intn, data_A, intm, slice_ptr_jpvt, slice_ptr_tau, work, lwork, &info);
            } else {
                lapack::geqp3(intm, intn, data_A, intm, slice_ptr_jpvt, slice_ptr_tau, work, lwork, rwork.get(), &info);
            }
            for (CBLAS_INT i = 0; i < intn; i++) {
                slice_ptr_jpvt[i] -= 1; // geqp3 returns a 1-based index array, so subtract 1
            }
        }

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            return BatchStatus::ok;
        }

        // `R` is always required, the correct shape is ensured by `middle_dim`.
        // NB. if a new array is allocated, it is pre-filled with zeros so only the relevant triangle should be copied.
        if (mode == QR_mode::R && overwrite_a) { // Data in-place, but need to remove the other part of the factorization still to obtain upper-triangular matrix
            zero_other_triangle('U', data_A, middle_dim, intn, intm);
        } else if (mode == QR_mode::FULL && overwrite_a) {
            // The `R` array should be copied back into the `A` array, which is F-ordered, so pretend the transpose got copied.
            // Then zero out the other triangle since the original input data is still in place.
            copy_triangle_to_C(slice_ptr_A, data_A, intn, intm, intm * (npy_intp)sizeof(T), (npy_intp)sizeof(T), 'L');
            zero_other_triangle('U', slice_ptr_A, intm, intn, intm);
        } else {
            copy_triangle_to_C(slice_ptr_R, data_A, middle_dim, intn, (npy_intp)sizeof(T), intm * (npy_intp)sizeof(T), 'U'); // F-ordered input to `copy_triangle`, hence `s0 == sizeof(T)` and `s1 == intm * sizeof(T)`
        }

        // Construct the Q matrix if required, same handling for both modes; dimensions are set already.
        if (mode == QR_mode::FULL || mode == QR_mode::ECONOMIC) {
            if constexpr (!is_complex_v<T>) {
                lapack::orgqr(intm, middle_dim, K, data_A, intm, slice_ptr_tau, work, lwork, &info);
            } else {
                lapack::ungqr(intm, middle_dim, K, data_A, intm, slice_ptr_tau, work, lwork, &info);
            }

            if (info != 0) {
                slice_status.lapack_info = (Py_ssize_t)info;
                vec_status.push_back(slice_status);
                return BatchStatus::ok;
            }

            if (!overwrite_a || mode == QR_mode::FULL) {
                copy_slice_F_to_C(slice_ptr_Q, data_A, intm, middle_dim);
            } // NB. else the data is already in place
        }
        // NB. for mode == RAW_MODE the data was either processed in place (for `overwrite_a`) or it was processed
        // directly in the return object. Hence, Q is already in place and no additional copy is needed.

    } // End of batching loop

    return BatchStatus::ok;
}

} // namespace sp_linalg
