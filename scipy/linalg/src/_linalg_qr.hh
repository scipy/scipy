/*
 * Templated loops for `linalg.qr`
 */

 #include "src/_common_array_utils.hh"
template<typename T>
int
_qr(PyArrayObject *ap_Am, PyArrayObject *ap_Q, PyArrayObject *ap_R, PyArrayObject *ap_tau, PyArrayObject *ap_jpvt, int overwrite_a, QR_mode mode, int pivoting, SliceStatusVec &vec_status)
{
    using real_type = typename type_traits<T>::real_type;
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

    if (mode == QR_mode::RAW) {
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
    CBLAS_INT middle_dim = (mode == QR_mode::ECONOMIC || mode == QR_mode::RAW) ? K : intm;


    // Probe both the factorization as well as `or_un_gqr` to find the optimal lwork
    T tmp_factor = numeric_limits<T>::zero;
    T tmp_or_un_gqr = numeric_limits<T>::zero;
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

    // M * max_dim for the buffer containing the data for A/Q.
    // If mode == RAW, the `tau` have to be returned, else a temporary buffer is sufficient.
    CBLAS_INT tau_size = (mode == QR_mode::RAW) ? 0 : K;
    T *buffer = (T *)malloc((M * max_dim + lwork + tau_size) * sizeof(T));
    if ( buffer == NULL ) { info = -101; return int(info); }

    T *data_A = &buffer[0];
    T *work = &buffer[M * max_dim];
    T *tau_buffer = &buffer[M * max_dim + lwork];

    // `c/zgeqp3` needs rwork
    void *rwork = NULL;
    if (pivoting && type_traits<T>::is_complex) {
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
        // Hence, the use of `ndim` instead of `ndim-1`.
        slice_ptr_A = compute_slice_ptr(idx, A_data, ndim, shape, strides);
        slice_ptr_R = compute_slice_ptr(idx, R_data, ndim, shape, strides_R);

        if (mode != QR_mode::R) {
            slice_ptr_Q = compute_slice_ptr(idx, Q_data, ndim, shape, strides_Q);
        }
        if (mode == QR_mode::RAW) {
            slice_ptr_tau = compute_slice_ptr(idx, tau_data, ndim, shape, strides_tau);
        } else { // `tau` might still be required, but should not be returned, so buffer is sufficient.
            slice_ptr_tau = tau_buffer;
        }
        if (pivoting) {
            slice_ptr_jpvt = compute_slice_ptr(idx, jpvt_data, ndim, shape, strides_jpvt);
        }

        copy_slice_F(data_A, slice_ptr_A, M, N, strides[ndim-2], strides[ndim-1]);

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

        // Extract useful information and continue computation if needed.
        // `R` is always required, the correct shape is ensured by `middle_dim`.
        extract_upper_triangle(slice_ptr_R, data_A, middle_dim, intn, intm);

        // Construct the Q matrix if required, same handling for both modes; dimensions are set already.
        if (mode == QR_mode::FULL || mode == QR_mode::ECONOMIC) {
            call_or_un_gqr(&intm, &middle_dim, &K, data_A, &intm, slice_ptr_tau, work, &lwork, &info);

            if (info != 0) {
                slice_status.lapack_info = (Py_ssize_t)info;
                vec_status.push_back(slice_status);
                goto free_exit;
            }

            copy_slice_F_to_C(slice_ptr_Q, data_A, intm, middle_dim);
        }

        // In the case of `raw` mode, the output of `geqrf`/`geqp3` is what is expected in `Q`.
        if (mode == QR_mode::RAW) {
            copy_slice_F_to_C(slice_ptr_Q, data_A, intm, intn);
        }

    } // End of batching loop

free_exit:
    free(buffer);
    free(rwork);

    return 1;
}
