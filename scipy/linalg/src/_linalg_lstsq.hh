/*
 * Templated loops for linalg.lstsq
 */

template<typename T>
int
_lstsq_gelss(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_S, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, SliceStatusVec& vec_status)
{
    using real_type = typename type_traits<T>::real_type; // float if T==npy_cfloat etc
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
    T tmp = numeric_limits<T>::zero;
    call_gelss(&intm, &intn, &int_nrhs, NULL, &lda, NULL, &ldb, NULL, &r_rcond, &rank, &tmp, &lwork, NULL, &info);
    if(info != 0) { return -100; }

    lwork = _calc_lwork(tmp);
    if (lwork < 0) {return -111;}

    // allocate
    npy_intp bufsize = m*n + ldb*nrhs + lwork;

    T* buffer = (T *)malloc((bufsize)*sizeof(T));
    if (NULL == buffer) { return -101; }

    // Chop the buffer into parts
    T* data = &buffer[0];
    T *data_b = &buffer[m*n];
    T *work = &buffer[m*n + ldb*nrhs];

    real_type *rwork = NULL;
    if constexpr (type_traits<T>::is_complex) {
        rwork = (real_type *)malloc(5*min_mn*sizeof(real_type));

        if (rwork == NULL) {
            free(buffer);
            return -102;
        }
    }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++){
        init_status(slice_status, idx, St::GENERAL);

        // copy the slice to `data` in F order
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
        copy_slice_F(data, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);

        // copy the r.h.s, too; NB: gelss needs LDB=max(1, m, n)
        T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim_b, shape_b, strides_b);
        copy_slice_F(data_b, slice_ptr_b, m, nrhs, strides_b[ndim_b-2], strides_b[ndim_b-1], ldb);

        // perform the least squares
        call_gelss(&intm, &intn, &int_nrhs, data, &lda, data_b, &ldb, ptr_S, &r_rcond, &rank, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy results from temp buffers (S is filled in-place already)
        copy_slice_F_to_C(ptr_x, data_b, n, nrhs, ldb);
        *ptr_rank = (npy_int64)rank;
        // NB: we discard the column residuals, b[n:]

        // advance the output pointers: S, x and rank arrays are C-ordered by construction
        ptr_S += min_mn;
        ptr_x += n*nrhs;
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
_lstsq_gelsd(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_S, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, SliceStatusVec& vec_status)
{
    using real_type = typename type_traits<T>::real_type; // float if T==npy_cfloat etc
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
    T tmp = numeric_limits<T>::zero;
    real_type tmp_lrwork = 0;
    CBLAS_INT liwork = 0, lrwork = 0;
    call_gelsd(&intm, &intn, &int_nrhs, NULL, &lda, NULL, &ldb, NULL, &r_rcond, &rank, &tmp, &lwork, &tmp_lrwork, &liwork, &info);

    if(info != 0) { return -100; }
    lwork = _calc_lwork(tmp);
    lrwork = type_traits<T>::is_complex ? _calc_lwork(tmp_lrwork) : 0 ; 

    if ((lwork < 0) || (lrwork < 0) || (liwork < 0)) {return -111;}

    // allocate
    npy_intp bufsize = m*n + ldb*nrhs + lwork;

    T* buffer = (T *)malloc((bufsize)*sizeof(T));
    if (NULL == buffer) { return -101; }

    // Chop the buffer into parts
    T* data = &buffer[0];
    T *data_b = &buffer[m*n];
    T *work = &buffer[m*n + ldb*nrhs];

    real_type *rwork = NULL;
    if constexpr (type_traits<T>::is_complex) {
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

        // copy the slice to `data` in F order
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
        copy_slice_F(data, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);

        // copy the r.h.s, too; NB: gelsd needs LDB=max(1, m, n)
        T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim_b, shape_b, strides_b);
        copy_slice_F(data_b, slice_ptr_b, m, nrhs, strides_b[ndim_b-2], strides_b[ndim_b-1], ldb);

        // perform the least squares
        call_gelsd(&intm, &intn, &int_nrhs, data, &lda, data_b, &ldb, ptr_S, &r_rcond, &rank, work, &lwork, rwork, iwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy results from temp buffers (S is filled in-place already)
        copy_slice_F_to_C(ptr_x, data_b, n, nrhs, ldb);
        *ptr_rank = (npy_int64)rank;
        // XXX we discard the column residuals, b[n:]

        // advance the output pointers: S, x and rank arrays are C-ordered by construction
        ptr_S += min_mn;
        ptr_x += n*nrhs;
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
_lstsq_gelsy(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, SliceStatusVec& vec_status)
{
    using real_type = typename type_traits<T>::real_type; // float if T==npy_cfloat etc
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
    T tmp = numeric_limits<T>::zero;
    call_gelsy(&intm, &intn, &int_nrhs, NULL, &lda, NULL, &ldb, NULL, &r_rcond, &rank, &tmp, &lwork, NULL, &info);
    if(info != 0) { return -100; }

    lwork = _calc_lwork(tmp);
    if (lwork < 0) {return -111;}

    // allocate
    npy_intp bufsize = m*n + ldb*nrhs + lwork;

    T* buffer = (T *)malloc((bufsize)*sizeof(T));
    if (NULL == buffer) { return -101; }

    // Chop the buffer into parts
    T* data = &buffer[0];
    T *data_b = &buffer[m*n];
    T *work = &buffer[m*n + ldb*nrhs];

    real_type *rwork = NULL;
    if constexpr (type_traits<T>::is_complex) {
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

        // copy the slice to `data` in F order
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
        copy_slice_F(data, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);

        // copy the r.h.s, too; NB: gelss needs LDB=max(1, m, n)
        T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim_b, shape_b, strides_b);
        copy_slice_F(data_b, slice_ptr_b, m, nrhs, strides_b[ndim_b-2], strides_b[ndim_b-1], ldb);

        for(npy_intp i=0; i<n; i++) {jpvt[i] = 0;}

        // perform the least squares
        call_gelsy(&intm, &intn, &int_nrhs, data, &lda, data_b, &ldb, jpvt, &r_rcond, &rank, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy results from temp buffers
        copy_slice_F_to_C(ptr_x, data_b, n, nrhs, ldb);
        *ptr_rank = (npy_int64)rank;
        // XXX we discard the column residuals, b[n:]

        // advance the output pointers: x and rank arrays are C-ordered by construction
        ptr_x += n*nrhs;
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
_lstsq(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_S, PyArrayObject *ap_x, PyArrayObject *ap_rank, double rcond, const char * lapack_driver, SliceStatusVec& vec_status)
{
    int info;
    if (strcmp(lapack_driver, "gelss") == 0) {
        info = _lstsq_gelss<T>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, vec_status);
    }
    else if (strcmp(lapack_driver, "gelsd") == 0) {
        info = _lstsq_gelsd<T>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, vec_status);
    }
    else if (strcmp(lapack_driver, "gelsy") == 0) {
        info = _lstsq_gelsy<T>(ap_Am, ap_b, ap_x, ap_rank, rcond, vec_status);
    }
    else {
        // should have been validated at call site, really
        info = -110;
    }
    return info;
}
