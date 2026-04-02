/*
 * Templated loops for `linalg.qr`
 */

template<typename T>
int
_qr(PyArrayObject *ap_Am, PyArrayObject *ap_Q, PyArrayObject *ap_R, PyArrayObject *ap_tau, PyArrayObject *ap_jpvt, int overwrite_a, QR_mode mode, int pivoting, SliceStatusVec &vec_status)
{
    using real_type = typename sp_type_traits<T>::real_type;
    SliceStatus slice_status;

    // ------------------------------------------------------------------------
    // Input array attributes, some bookkeeping due to conditional allocations
    // ------------------------------------------------------------------------
    T *A_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp M = shape[ndim-2];
    npy_intp N = shape[ndim-1];

    T *R_data = (T *)PyArray_DATA(ap_R);
    npy_intp *strides_R = PyArray_STRIDES(ap_R);

    // Conditional allocations based on the mode; not all data is always available.
    T *Q_data = NULL, *tau_data = NULL;
    CBLAS_INT *jpvt_data = NULL;
    npy_intp *strides_Q = NULL, *strides_tau = NULL, *strides_jpvt = NULL;

    if (mode != QR_mode::R) {
        Q_data = (T *)PyArray_DATA(ap_Q);
        strides_Q = PyArray_STRIDES(ap_Q);
    }

    if (mode == QR_mode::RAW_MODE) {
        tau_data = (T *)PyArray_DATA(ap_tau);
        strides_tau = PyArray_STRIDES(ap_tau);
    }

    if (pivoting) {
        jpvt_data = (CBLAS_INT *)PyArray_DATA(ap_jpvt);
        strides_jpvt = PyArray_STRIDES(ap_jpvt);
    }

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    // -------------------------------------------------------------------
    // Workspace computation and allocation
    // -------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)N, intm = (CBLAS_INT)M, info = 0;
    CBLAS_INT K = std::min(intn, intm), max_dim = std::max(intn, intm);
    CBLAS_INT middle_dim = (mode == QR_mode::ECONOMIC || mode == QR_mode::RAW_MODE) ? K : intm; // Final dimension: Q is [`M` x `middle_dim`], R is [`middle_dim` x `N`]


    // Probe both the factorization as well as `or_un_gqr` to find the optimal lwork
    T tmp_factor = sp_numeric_limits<T>::zero;
    T tmp_or_un_gqr = sp_numeric_limits<T>::zero;
    CBLAS_INT lwork = -1;

    if (!pivoting) {
        call_geqrf(&intm, &intn, NULL, &intm, NULL, &tmp_factor, &lwork, &info);
    } else {
        call_geqp3(&intm, &intn, NULL, &intm, NULL, NULL, &tmp_factor, &lwork, NULL, &info);
    }
    if (info != 0) { info = -100; return (int)info; }

    // Also probe the Q construction step if required.
    if (mode == QR_mode::FULL || mode == QR_mode::ECONOMIC) {
        call_or_un_gqr(&intm, &middle_dim, &K, NULL, &intm, NULL, &tmp_or_un_gqr, &lwork, &info);
    }
    if (info != 0) { info = -100; return (int)info; }

    CBLAS_INT lwork_factor = _calc_lwork(tmp_factor);
    CBLAS_INT lwork_or_un_gqr = _calc_lwork(tmp_or_un_gqr);
    lwork = std::max(lwork_factor, lwork_or_un_gqr);

    if ((lwork < 0) || (pivoting && 3 * N + 1 >= std::numeric_limits<CBLAS_INT>::max())) {
        // Too large lwork required - Computation cannot be performed.
        // `geqrf` and `or_un_gqr` require an `lwork` of at least N,
        // `geqp3` requires an `lwork` of at least 3N + 1. The former would not
        // be detectable, but the latter is.
        return -99;
    }

    /*
     * Allocate buffer and slice into parts
     *
     *     lwork        tau_size     buf_size
     * |------------|-------------|--------------|
     * ^            ^             ^
     * work         tau_buffer    data_A
     *
     * - `work` is the work area, size `lwork`
     * - `tau_buffer` contains a temporary buffer for the reflectors. If mode == RAW_MODE, they have
     *   to be returned and hence a buffer is not needed. If needed, size = K
     * - `data_A` contains a temporary buffer for LAPACK calls, if mode == FULL, the resulting
     *   `Q` is `M` x `M`. However, `max_dim` is used to allocate enough space should `M` < `N`.
     *   If `overwrite_a` is enabled this buffer is not needed. Similarly, if mode == RAW the
     *   returned argument is also in F-order, so the intermediate buffer step can be skipped.
     */
    CBLAS_INT tau_size = (mode == QR_mode::RAW_MODE) ? 0 : K;
    CBLAS_INT buf_size = ((overwrite_a && mode != QR_mode::FULL) || mode == QR_mode::RAW_MODE) ? 0 : (mode == QR_mode::FULL) ? M * max_dim : M * N;
    T *buffer = (T *)malloc((buf_size + lwork + tau_size) * sizeof(T));
    if ( buffer == NULL ) { info = -101; return int(info); }

    T *work = &buffer[0];
    T *tau_buffer = &buffer[lwork];
    T *data_A = &buffer[lwork + tau_size];

    // `c/zgeqp3` needs rwork
    void *rwork = NULL;
    if (pivoting && sp_type_traits<T>::is_complex) {
        rwork = malloc(2 * N * sizeof(real_type));

        if (rwork == NULL) {
            free(buffer);
            info = -102;
            return int(info);
        }
    }

    // -------------------------------------------------------------------
    // Main loop to traverse the slices
    // -------------------------------------------------------------------
    T *slice_ptr_A = NULL, *slice_ptr_Q = NULL, *slice_ptr_R = NULL, *slice_ptr_tau = NULL;
    CBLAS_INT *slice_ptr_jpvt = NULL;
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        // Compute slice pointers if necessary. `tau` and `jpvt` are vectors, so actually 1 dimension
        // less, however, the looping to compute the slices should traverse the same number of dimensions.
        // Hence, the use of `ndim` instead of `ndim-1`. Similarly, the shapes of `Q` and `R` might differ,
        // but only the batching ones are relevant, which are identical.
        slice_ptr_A = compute_slice_ptr(idx, A_data, ndim, shape, strides);
        if (!(mode == QR_mode::R && overwrite_a)) {
            slice_ptr_R = compute_slice_ptr(idx, R_data, ndim, shape, strides_R);
        } // NB. if the condition is true the relevant buffer is already in place

        if (mode != QR_mode::R) {
            slice_ptr_Q = compute_slice_ptr(idx, Q_data, ndim, shape, strides_Q);
        }
        if (mode == QR_mode::RAW_MODE) {
            slice_ptr_tau = compute_slice_ptr(idx, tau_data, ndim, shape, strides_tau);
        } else { // `tau` might still be required, but should not be returned, so buffer is sufficient.
            slice_ptr_tau = tau_buffer;
        }
        if (pivoting) {
            slice_ptr_jpvt = compute_slice_ptr(idx, jpvt_data, ndim, shape, strides_jpvt);
        }

        /*
         * Copy the buffer for processing. When `overwrite_a` is disabled the data should be handled in
         * another buffer, so copying is necessary. When it is set, the data can be processed in place.
         *
         * Another optimization is to process the data for mode == RAW in the return object directly as
         * the returned array is in F-order, hence bypassing the buffer avoids two redundant copies.
         */
        if ((!overwrite_a && mode != QR_mode::RAW_MODE) || mode == QR_mode::FULL) {
            copy_slice_F(data_A, slice_ptr_A, M, N, strides[ndim-2], strides[ndim-1]);
        } else if (!overwrite_a && mode == QR_mode::RAW_MODE) {
            copy_slice_F(slice_ptr_Q, slice_ptr_A, M, N, strides[ndim-2], strides[ndim-1]);
            data_A = slice_ptr_Q; // Add alias for easier to read codepaths later on.
        } else if (overwrite_a && (mode == QR_mode::RAW_MODE || mode == QR_mode::R || mode == QR_mode::ECONOMIC)) {
            data_A = slice_ptr_A;
        }

        init_status(slice_status, idx, St::GENERAL);

        // ---------------------------------------------------------------
        // Heavy lifting: call appropriate LAPACK routines.
        // ---------------------------------------------------------------

        // Factorization step, identical for each of the modes so do jointly.
        if (pivoting) {
            call_geqp3(&intm, &intn, data_A, &intm, slice_ptr_jpvt, slice_ptr_tau, work, &lwork, rwork, &info);
            for (CBLAS_INT i = 0; i < intn; i++) {
                slice_ptr_jpvt[i] -= 1; // geqp3 returns a 1-based index array, so subtract 1
            }
        } else {
            call_geqrf(&intm, &intn, data_A, &intm, slice_ptr_tau, work, &lwork, &info);
        }

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            goto free_exit;
        }

        // `R` is always required, the correct shape is ensured by `middle_dim`.
        if (mode == QR_mode::R && overwrite_a) {
            zero_other_triangle('U', data_A, middle_dim, intn, intm); // Data in-place, but need to remove the other part of the factorization still
        } else if (mode == QR_mode::FULL && overwrite_a) {
            // The `R` array should be copied back into the `A` array, which is F-ordered, so pretend the transpose got copied.
            // Then zero out the other triangle since the original input data is still in place.
            copy_triangle_to_C(slice_ptr_A, data_A, intn, intm, intm, 1, 'L');
            zero_other_triangle('U', slice_ptr_A, intm, intn, intm);
        } else {
            copy_triangle_to_C(slice_ptr_R, data_A, middle_dim, intn, 1, intm, 'U'); // F-ordered input to `copy_triangle`, hence `s0 == 1` and `s1 == intm`
        } // NB. if a new array is allocated, it is pre-filled with zeros so only the relevant triangle should be copied.

        // Construct the Q matrix if required, same handling for both modes; dimensions are set already.
        if (mode == QR_mode::FULL || mode == QR_mode::ECONOMIC) {
            call_or_un_gqr(&intm, &middle_dim, &K, data_A, &intm, slice_ptr_tau, work, &lwork, &info);

            if (info != 0) {
                slice_status.lapack_info = (Py_ssize_t)info;
                vec_status.push_back(slice_status);
                goto free_exit;
            }

            if (!overwrite_a || mode == QR_mode::FULL) {
                copy_slice_F_to_C(slice_ptr_Q, data_A, intm, middle_dim);
            } // NB. else the data is already in place
        }
        // NB. for mode == RAW the data was either processed in place (for `overwrite_a`),
        // or was processed in the return object. Hence, Q is already in place and no additional copy is needed.

    } // End of batching loop

free_exit:
    free(buffer);
    free(rwork);

    return 1;
}
