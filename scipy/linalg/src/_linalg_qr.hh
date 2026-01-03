/*
 * Templated loops for `linalg.qr`
 */
 #include "Python.h"
 #include <iostream>
 #include "numpy/arrayobject.h"
 #include "numpy/npy_math.h"
 #include "npy_cblas.h"
 #include "_npymath.hh"
 #include "_common_array_utils.hh"

template<typename T>
int
_qr(PyArrayObject *ap_Am, PyArrayObject *ap_Q, PyArrayObject *ap_R, PyArrayObject *ap_P, int overwrite_a, int lwork, QR_mode mode, int pivot, SliceStatusVec &vec_status)
{
    using real_type = typename type_traits<T>::real_type;

    SliceStatus slice_status;

    // -------------------------------------------------------------------
    // Input array attributes
    // -------------------------------------------------------------------
    T *A_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp M = shape[ndim-2];
    npy_intp N = shape[ndim-1];

    T *Q_data = (T *)PyArray_DATA(ap_Q);
    npy_intp *strides_Q = PyArray_STRIDES(ap_Q);
    T *R_data = (T *)PyArray_DATA(ap_R);
    npy_intp *strides_R = PyArray_STRIDES(ap_R);
    CBLAS_INT *P_data = (CBLAS_INT *)PyArray_DATA(ap_P);
    npy_intp *strides_P = PyArray_STRIDES(ap_P);

    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }
    }

    // -------------------------------------------------------------------
    // Workspace computation and allocation
    // -------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)N, intm = (CBLAS_INT)M, info = 0;
    CBLAS_INT K = std::min(intn, intm);

    // If the size of the array is not given, first do a "probe" call to
    // compute it. Logic is taken from `_linalg_solve()`.
    if (lwork == -1) {
        T tmp = numeric_limits<T>::zero;
        geqrf(&intm, &intn, NULL, &intm, NULL, &tmp, &lwork, &info);
        if (info != 0) { info = -100; return (int)info; }

        lwork = _calc_lwork(tmp);
        if ((lwork < 0) ||
            (N > std::numeric_limits<CBLAS_INT>::max() / 3)
        ) {
            // Too large lwork required - Computation cannot be performed.
            return -99;
        }

        // Guards for the calls to LAPACK routines
        lwork = (3 * N + 1 > lwork ? 3 * N + 1 : lwork);  // `s/dgeqp3` needs lwork of at least 3n + 1
        lwork = (M > lwork ? M : lwork);  // guard for the call to `or/ungqr`.
    }
    else if (!pivot && lwork < std::max(M, N)) {
        // The assigned `lwork` would not be enough to perform computations
        PyErr_SetString(PyExc_ValueError, "Without pivoting an lwork of at least max(M, N) is required.");
        return -1;
    }
    else if (pivot && lwork < std::max(M, 3 * N + 1)) {
        PyErr_SetString(PyExc_ValueError, "With pivoting and real arrays, an lwork of at least max(M, 3 * N + 1) is required.");
        return -1;
    }


    // `std::min(M, N)` for `tau`, std::max(M, N) for the buffer for A/Q
    T *buffer = (T *)malloc((M * std::max(M, N) + lwork + std::min(M, N)) * sizeof(T));
    if ( buffer == NULL ) { info = -101; return int(info); }

    T *data_A = &buffer[0];
    T *work = &buffer[M * std::max(M, N)];
    T *tau = &buffer[M * std::max(M, N) + lwork];

    // `c/zgeqp3` needs rwork
    void *rwork = NULL;
    if (type_traits<T>::is_complex) {
        rwork = malloc(2 * N * sizeof(real_type));

        if (rwork == NULL) {
            free(buffer);
            info = -102;
            return int(info);
        }
    }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        // Bundle all looping for the slice pointers into one large loop.
        // The shape of all matrices is the same across the batching dimensions.
        npy_intp offset_A = 0, offset_Q = 0, offset_R = 0, offset_P = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset_A += (temp_idx % shape[i]) * strides[i];
            offset_Q += (temp_idx % shape[i]) * strides_Q[i];
            offset_R += (temp_idx % shape[i]) * strides_R[i];
            offset_P += (temp_idx % shape[i]) * strides_P[i];

            temp_idx /= shape[i];
        }

        T *slice_ptr_A = (T *)(A_data + (offset_A / sizeof(T)));
        T *slice_ptr_Q = (T *)(Q_data + (offset_Q / sizeof(T)));
        T *slice_ptr_R = (T *)(R_data + (offset_R / sizeof(T)));
        CBLAS_INT *slice_ptr_P = (CBLAS_INT *)(P_data + (offset_P/sizeof(CBLAS_INT)));

        copy_slice_F(data_A, slice_ptr_A, M, N, strides[ndim-2], strides[ndim-1]);

        init_status(slice_status, idx, St::GENERAL);

        // Factorization step is identical for each algorithm so do not delegate
        if (pivot) {
            memset(slice_ptr_P, 0, intn * sizeof(int)); // geqp3 also takes in pivoting elements
            geqp3(&intm, &intn, data_A, &intm, slice_ptr_P, tau, work, &lwork, rwork, &info);
            for (int i = 0; i < intn; i++) {
                slice_ptr_P[i] -= 1; // geqp3 returns a 1-based index array, so subtract 1
            }
        }
        else {
            geqrf(&intm, &intn, data_A, &intm, tau, work, &lwork, &info);
        }

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            goto free_exit;
        }

        // TODO: accomodate for different solutions
        switch (mode) {
            case QR_mode::FULL:
            {
                extract_upper_triangle(slice_ptr_R, data_A, intm, intn, intm);

                // Full mode QR, hence Q will be MxM.
                // N.B. the number of reflectors is limited by the smallest dimension (= `K`)
                orungqr(&intm, &intm, &K, data_A, &intm, tau, work, &lwork, &info);

                if (info != 0) {
                    slice_status.lapack_info = (Py_ssize_t)info;
                    vec_status.push_back(slice_status);

                    goto free_exit;
                }

                copy_slice_F_to_C(slice_ptr_Q, data_A, intm, intm);
                break;
            }

            case QR_mode::R:
            {
                break;
            }

            case QR_mode::RAW:
            {
                break;
            }

            case QR_mode::ECONOMIC:
            {
                extract_upper_triangle(slice_ptr_R, data_A, K, intn, intm);

                // Economic mode QR, hence Q is MxK
                orungqr(&intm, &K, &K, data_A, &intm, tau, work, &lwork, &info);

                if (info != 0) {
                    slice_status.lapack_info = (Py_ssize_t)info;
                    vec_status.push_back(slice_status);

                    goto free_exit;
                }

                copy_slice_F_to_C(slice_ptr_Q, data_A, intm, K);
                break;
            }
        }
    }

free_exit:
    free(buffer);
    free(rwork);

    return 1;
}
