/*
 * Templated loops for linalg.lstsq
 */

namespace sp_linalg {

template<typename T>
int
_lstsq_gelss(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_S, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, const int overwrite_a, const int overwrite_b, SliceStatusVec& vec_status)
{
    using real_type = typename detail::type_traits<T>::real_type; // float if T==npy_cfloat etc
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
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    T *bm_data = (T *)PyArray_DATA(ap_b);
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);
    npy_intp *strides_b = PyArray_STRIDES(ap_b);
    npy_intp nrhs = PyArray_DIM(ap_b, ndim_b -1); // Number of right-hand-sides

    // XXX: the call site ensured that
    // 1. a.ndim == b.ndim
    // 2. a.shape = batch_shape + (m, n), and
    //    b.shape = batch_shape + (m, nrhs)
    // so we are not checking this again.

    // Outputs
    real_type *ptr_S = (real_type *)PyArray_DATA(ap_S);
    T *ptr_x = (T *)PyArray_DATA(ap_x);
    npy_int64 *ptr_rank = (npy_int64 *)PyArray_DATA(ap_rank);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, intm = (CBLAS_INT)m, int_nrhs = (CBLAS_INT)nrhs;
    CBLAS_INT lwork = -1, info;
    CBLAS_INT lda = intm > 1 ? intm : 1;

    CBLAS_INT max_mn = intm > intn ? intm : intn;
    CBLAS_INT min_mn = intm > intn ? intn : intm;
    CBLAS_INT ldb = max_mn > 1 ? max_mn : 1;
    real_type r_rcond = (real_type)rcond;
    CBLAS_INT rank = min_mn;

    // query LWORK
    T tmp = detail::numeric_limits<T>::zero;
    call_gelss(&intm, &intn, &int_nrhs, NULL, &lda, NULL, &ldb, NULL, &r_rcond, &rank, &tmp, &lwork, NULL, &info);
    if(info != 0) { return -100; }

    lwork = _calc_lwork(tmp);
    if (lwork < 0) {return -111;}

    /*
     * Allocate buffer and chop up into parts
     *
     *   lwork     data_a_size   data_b_size
     * |---------|-------------|-------------|
     * ^         ^             ^
     * work      data_a        data_b
     *
     * - `work` is the work area, size `lwork`
     * - `data_a` is the buffer containing the `a` matrix, size `m * n` if `overwrite_a` disabled, else 0.
     * - `data_b` is the buffer containing the `b` matrix, size `ldb * nrhs` else 0.
     */
    npy_intp data_a_size = overwrite_a ? 0 : m * n;
    npy_intp data_b_size = overwrite_b ? 0 : ldb * nrhs;
    npy_intp bufsize = data_a_size + data_b_size + lwork;

    T* buffer = (T *)malloc((bufsize)*sizeof(T));
    if (NULL == buffer) { return -101; }

    // Chop the buffer into parts
    // NB. `overwrite_x` only enabled for 2D F-contiguous arrays
    T *work = &buffer[0];

    T *data_a = NULL;
    if (!overwrite_a) {
        data_a = &buffer[lwork];
    } else {
        data_a = Am_data;
    }

    T *data_b = NULL;
    if (!overwrite_b) {
        data_b = &buffer[lwork + data_a_size];
    } else {
        data_b = bm_data;
    }

    real_type *rwork = NULL;
    if constexpr (detail::type_traits<T>::is_complex) {
        rwork = (real_type *)malloc(5*min_mn*sizeof(real_type));

        if (rwork == NULL) {
            free(buffer);
            return -102;
        }
    }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++){
        init_status(slice_status, idx, St::GENERAL);

        // copy the slices to `data` in F order. If `overwrite_x` is enabled the data
        // is already in place, so no need to.
        if (!overwrite_a) {
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_slice_F(data_a, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);
        }

        if (!overwrite_b) {
            T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim_b, shape_b, strides_b);
            copy_slice_F(data_b, slice_ptr_b, m, nrhs, strides_b[ndim_b-2], strides_b[ndim_b-1], ldb);
        } // NB. gelss needs LDB = max(1, m, n)

        // perform the least squares
        call_gelss(&intm, &intn, &int_nrhs, data_a, &lda, data_b, &ldb, ptr_S, &r_rcond, &rank, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        if (!overwrite_b) {
            // copy results from temp buffers (S is filled in-place already)
            copy_slice_F_to_C(ptr_x, data_b, ldb, nrhs, ldb);
        } // NB. if `overwrite_b` is true the data is already in place and should just be shrunk down at the python side.
        *ptr_rank = (npy_int64)rank;

        // advance the output pointers: S, x and rank arrays are C-ordered by construction
        ptr_S += min_mn;
        ptr_x += ldb*nrhs;
        ptr_rank += 1;
    }

done:
    free(buffer);
    free(rwork);

    return 1;
}


// XXX consider deduplicating _gelsd, gelss and gelss loops.

template<typename T>
int
_lstsq_gelsd(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_S, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, const int overwrite_a, const int overwrite_b, SliceStatusVec& vec_status)
{
    using real_type = typename detail::type_traits<T>::real_type; // float if T==npy_cfloat etc
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
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    T *bm_data = (T *)PyArray_DATA(ap_b);
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);
    npy_intp *strides_b = PyArray_STRIDES(ap_b);
    npy_intp nrhs = PyArray_DIM(ap_b, ndim_b -1); // Number of right-hand-sides

    // Outputs
    real_type *ptr_S = (real_type *)PyArray_DATA(ap_S);
    T *ptr_x = (T *)PyArray_DATA(ap_x);
    npy_int64 *ptr_rank = (npy_int64 *)PyArray_DATA(ap_rank);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, intm = (CBLAS_INT)m, int_nrhs = (CBLAS_INT)nrhs;
    CBLAS_INT lwork = -1, info;
    CBLAS_INT lda = intm > 1 ? intm : 1;

    CBLAS_INT max_mn = intm > intn ? intm : intn;
    CBLAS_INT min_mn = intm > intn ? intn : intm;
    CBLAS_INT ldb = max_mn > 1 ? max_mn : 1;
    real_type r_rcond = (real_type)rcond;
    CBLAS_INT rank = min_mn;

    // query LWORK, LRWORK and LIWORK
    // XXX: bump LRWORK and LIWORK up to improve perf? lwork=-1 query returns the *minimum* values
    T tmp = detail::numeric_limits<T>::zero;
    real_type tmp_lrwork = 0;
    CBLAS_INT liwork = 0, lrwork = 0;
    call_gelsd(&intm, &intn, &int_nrhs, NULL, &lda, NULL, &ldb, NULL, &r_rcond, &rank, &tmp, &lwork, &tmp_lrwork, &liwork, &info);

    if(info != 0) { return -100; }
    lwork = _calc_lwork(tmp);
    lrwork = detail::type_traits<T>::is_complex ? _calc_lwork(tmp_lrwork) : 0 ;

    if ((lwork < 0) || (lrwork < 0) || (liwork < 0)) {return -111;}

    /*
     * Allocate buffer and chop up into parts
     *
     *   lwork     data_a_size   data_b_size
     * |---------|-------------|-------------|
     * ^         ^             ^
     * work      data_a        data_b
     *
     * - `work` is the work area, size `lwork`
     * - `data_a` is the buffer containing the `a` matrix, size `m * n` if `overwrite_a` disabled, else 0.
     * - `data_b` is the buffer containing the `b` matrix, size `ldb * nrhs` else 0.
     */
    npy_intp data_a_size = overwrite_a ? 0 : m * n;
    npy_intp data_b_size = overwrite_b ? 0 : ldb * nrhs;
    npy_intp bufsize = data_a_size + data_b_size + lwork;

    T* buffer = (T *)malloc((bufsize)*sizeof(T));
    if (NULL == buffer) { return -101; }

    // Chop the buffer into parts
    // NB. `overwrite_x` is only enabled for 2D F-contiguous arrays
    T *work = &buffer[0];

    T *data_a = NULL;
    if (!overwrite_a) {
        data_a = &buffer[lwork];
    } else {
        data_a = Am_data;
    }

    T *data_b = NULL;
    if (!overwrite_b) {
        data_b = &buffer[lwork + data_a_size];
    } else {
        data_b = bm_data;
    }

    real_type *rwork = NULL;
    if constexpr (detail::type_traits<T>::is_complex) {
        rwork = (real_type *)malloc(lrwork*sizeof(real_type));

        if (rwork == NULL) {
            free(buffer);
            return -102;
        }
    }

    CBLAS_INT *iwork = (CBLAS_INT *)malloc(liwork*sizeof(CBLAS_INT));
    if (iwork == NULL) {
        free(buffer);
        free(rwork);
        return -103;
    }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++){
        init_status(slice_status, idx, St::GENERAL);

        // copy data to buffers in F order. NB. `overwrite_x` is only enabled for 2D
        // F-contiguous data, so the input would already be in place.
        if (!overwrite_a) {
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_slice_F(data_a, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);
        }

        if (!overwrite_b) {
            // copy the r.h.s, too; NB: gelsd needs LDB=max(1, m, n)
            T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim_b, shape_b, strides_b);
            copy_slice_F(data_b, slice_ptr_b, m, nrhs, strides_b[ndim_b-2], strides_b[ndim_b-1], ldb);
        } // NB. gelsd needs ldb = max(1, m, n)

        // perform the least squares
        call_gelsd(&intm, &intn, &int_nrhs, data_a, &lda, data_b, &ldb, ptr_S, &r_rcond, &rank, work, &lwork, rwork, iwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        if (!overwrite_b) {
            // copy results from temp buffers (S is filled in-place already)
            copy_slice_F_to_C(ptr_x, data_b, ldb, nrhs, ldb);
        } // NB. else the data is already in place
        *ptr_rank = (npy_int64)rank;

        // advance the output pointers: S, x, residuals and rank arrays are C-ordered by construction
        ptr_S += min_mn;
        ptr_x += ldb*nrhs;
        ptr_rank += 1;
    }

done:
    free(buffer);
    free(rwork);
    free(iwork);
    return 1;
}


template<typename T>
int
_lstsq_gelsy(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, const int overwrite_a, const int overwrite_b, SliceStatusVec& vec_status)
{
    using real_type = typename detail::type_traits<T>::real_type; // float if T==npy_cfloat etc
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
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    T *bm_data = (T *)PyArray_DATA(ap_b);
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);
    npy_intp *strides_b = PyArray_STRIDES(ap_b);
    npy_intp nrhs = PyArray_DIM(ap_b, ndim_b -1); // Number of right-hand-sides

    // Outputs
    T *ptr_x = (T *)PyArray_DATA(ap_x);
    npy_int64 *ptr_rank = (npy_int64 *)PyArray_DATA(ap_rank);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, intm = (CBLAS_INT)m, int_nrhs = (CBLAS_INT)nrhs;
    CBLAS_INT lwork = -1, info;
    CBLAS_INT lda = intm > 1 ? intm : 1;

    CBLAS_INT max_mn = intm > intn ? intm : intn;
    CBLAS_INT min_mn = intm > intn ? intn : intm;
    CBLAS_INT ldb = max_mn > 1 ? max_mn : 1;
    real_type r_rcond = (real_type)rcond;
    CBLAS_INT rank = min_mn;

    // query LWORK
    T tmp = detail::numeric_limits<T>::zero;
    call_gelsy(&intm, &intn, &int_nrhs, NULL, &lda, NULL, &ldb, NULL, &r_rcond, &rank, &tmp, &lwork, NULL, &info);
    if(info != 0) { return -100; }

    lwork = _calc_lwork(tmp);
    if (lwork < 0) {return -111;}

    /*
     * Allocate buffer and chop up into parts
     *
     *   lwork     data_a_size   data_b_size
     * |---------|-------------|-------------|
     * ^         ^             ^
     * work      data_a        data_b
     *
     * - `work` is the work area, size `lwork`
     * - `data_a` is the buffer containing the `a` matrix, size `m * n` if `overwrite_a` disabled, else 0.
     * - `data_b` is the buffer containing the `b` matrix, size `ldb * nrhs` else 0.
     */
    npy_intp data_a_size = overwrite_a ? 0 : m * n;
    npy_intp data_b_size = overwrite_b ? 0 : ldb * nrhs;
    npy_intp bufsize = data_a_size + data_b_size + lwork;

    T* buffer = (T *)malloc((bufsize)*sizeof(T));
    if (NULL == buffer) { return -101; }

    T *work = &buffer[0];

    // NB. `overwrite_x` only enabled for 2D F-contiguous inputs
    T *data_a = NULL;
    if (!overwrite_a) {
        data_a = &buffer[lwork];
    } else {
        data_a = Am_data;
    }

    T *data_b = NULL;
    if (!overwrite_b) {
        data_b = &buffer[lwork + data_a_size];
    } else {
        data_b = bm_data;
    }

    real_type *rwork = NULL;
    if constexpr (detail::type_traits<T>::is_complex) {
        rwork = (real_type *)malloc(2*n*sizeof(real_type));

        if (rwork == NULL) {
            free(buffer);
            return -102;
        }
    }

    CBLAS_INT *jpvt = (CBLAS_INT *)malloc(n*sizeof(CBLAS_INT));
    if (jpvt == NULL) {
        free(buffer);
        free(rwork);
        return -103;
    }
    for(npy_intp i=0; i<n; i++) {jpvt[i] = 0;}


    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++){
        init_status(slice_status, idx, St::GENERAL);

        // copy the slices to buffers in F order. `overwrite_x` only enabled for 2D F-contiguous
        // inputs, so the data is already in place then.
        if (!overwrite_a) {
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_slice_F(data_a, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);
        }

        if (!overwrite_b) {
            T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim_b, shape_b, strides_b);
            copy_slice_F(data_b, slice_ptr_b, m, nrhs, strides_b[ndim_b-2], strides_b[ndim_b-1], ldb);
        } // NB. gelsy needs ldb = max(1, m, n)

        for(npy_intp i=0; i<n; i++) {jpvt[i] = 0;}

        // perform the least squares
        call_gelsy(&intm, &intn, &int_nrhs, data_a, &lda, data_b, &ldb, jpvt, &r_rcond, &rank, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        if (!overwrite_b) {
            // copy results from temp buffers
            copy_slice_F_to_C(ptr_x, data_b, ldb, nrhs, ldb);
        } // NB. Else the data is already in place
        *ptr_rank = (npy_int64)rank;

        // advance the output pointers: x and rank arrays are C-ordered by construction
        ptr_x += ldb*nrhs;
        ptr_rank += 1;
    }

done:
    free(buffer);
    free(rwork);
    free(jpvt);
    return 1;
}


template<typename T>
int
_lstsq(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_S, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, const char * lapack_driver, const int overwrite_a, const int overwrite_b, SliceStatusVec& vec_status)
{
    int info;
    if (strcmp(lapack_driver, "gelss") == 0) {
        info = _lstsq_gelss<T>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, overwrite_a, overwrite_b, vec_status);
    }
    else if (strcmp(lapack_driver, "gelsd") == 0) {
        info = _lstsq_gelsd<T>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, overwrite_a, overwrite_b, vec_status);
    }
    else if (strcmp(lapack_driver, "gelsy") == 0) {
        info = _lstsq_gelsy<T>(ap_Am, ap_b, ap_x, ap_rank, rcond, overwrite_a, overwrite_b, vec_status);
    }
    else {
        // should have been validated at call site, really
        info = -110;
    }
    return info;
}

} // namespace sp_linalg
