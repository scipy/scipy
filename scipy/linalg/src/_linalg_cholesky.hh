template<typename T>
int
_cholesky(PyArrayObject *ap_Am, PyArrayObject *ap_Cm, int lower, int overwrite_a, int clean, SliceStatusVec &vec_status) {
    char uplo = lower ? 'L' : 'U';
    St slice_structure = St::NONE;
    SliceStatus slice_status;

    //-------------------------------------------------------------------
    // Input array attributes
    //-------------------------------------------------------------------
    T *Am_data = (T *)PyArray_DATA(ap_Am);
    npy_intp ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp n = shape[ndim-1];

    T *ret_data = (T *)PyArray_DATA(ap_Cm); // Identical shape as input, C-ordered

    // Determine number of slices (np.prod(dim[:-2]))
    npy_intp outer_size = 1;
    if (ndim > 2) {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }
    }

    //--------------------------------------------------------------------
    // Workspace computation and allocation
    //--------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, info = 0;

    /*
     * Only need a buffer for `a`, but if `overwrite_a` is enabled this is not the case.
     *
     * No `work` buffer required for `potrf`.
     */
    npy_intp buf_size = overwrite_a ? 0 : n * n;
    T *buffer = (T *)malloc(buf_size * sizeof(T));
    if (buffer == NULL) { info = -104; return info; }

    T *data_a = NULL;
    if (!overwrite_a) {
        data_a = &buffer[0];
    } else {
        data_a = Am_data;
    }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        if (!overwrite_a) {
            T *slice_ptr_A = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_slice_F(data_a, slice_ptr_A, n, n, strides[ndim-2], strides[ndim-1]);
        }
        // NB. `overwrite_a` is only enabled when the input is 2D F-ordered so data is
        // already correctly ordered in input array. If generalized to `ndim > 2` this
        // should be updated accordingly.

        init_status(slice_status, idx, slice_structure);
        call_potrf(&uplo, &intn, data_a, &intn, &info);

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            goto done;
        }

        /*
         * Only copy over relevant part of the triangle. Other part is already
         * put to 0 by construction. For `overwrite_a` the story is a biff different
         * since the input array should then still be cleaned after the fact.
         *
         * NB. this `zero_other_triangle` only works since `overwrite_a` is only enabled
         * for 2D F-ordered arrays. Otherwise some other measures have to be taken.
         */
        if (!overwrite_a) {
            copy_triangle_F_to_C(&ret_data[idx * n * n], data_a, n, uplo);
        } else if (clean) {
            zero_other_triangle(uplo, data_a, n);
        }

    } // end of batching loop

done:
    free(buffer);
    return 1;
}
