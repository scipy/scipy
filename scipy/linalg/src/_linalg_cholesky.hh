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

    T *ret_data = (T *)PyArray_DATA(ap_Cm); // == Am_data when overwrite_a, else zero-initialized C-ordered

    // Determine number of slices (np.prod(dim[:-2]))
    npy_intp outer_size = 1;
    if (ndim > 2) {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }
    }

    //--------------------------------------------------------------------
    // Workspace computation and allocation
    //--------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, info = 0;

    npy_intp sz = (npy_intp)sizeof(T);
    bool f_contig = (strides[ndim-2] == sz && strides[ndim-1] == n * sz);

    /*
     * Since scipy returns C-ordered output, we call potrf with flipped uplo
     * so the result lands in the correct C-order triangle directly.
     *
     * For Hermitian input, the implicit conjugation from reinterpreting
     * C-order as F-order cancels with the conjugation from the uplo flip.
     *
     * overwrite_a with F-contiguous: potrf uses original uplo in-place.
     * overwrite_a with C-contiguous: potrf uses flipped uplo in-place.
     */
    npy_intp s0 = strides[ndim-2] / sz;
    npy_intp s1 = strides[ndim-1] / sz;

    // Flip uplo for all paths except F-contiguous overwrite_a.
    if (!(overwrite_a && f_contig)) {
        uplo = (uplo == 'L') ? 'U' : 'L';
    }

    T *data_a = overwrite_a ? Am_data : NULL;

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        if (!overwrite_a) {
            data_a = &ret_data[idx * n * n];
            T *slice_ptr_A = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_triangle_to_C(data_a, slice_ptr_A, n, s0, s1, uplo);
        }
        // NB. `overwrite_a` is only enabled when the input is 2D contiguous so
        // data is already correctly ordered in input array. If generalized to
        // `ndim > 2` this should be updated accordingly.

        init_status(slice_status, idx, slice_structure);
        call_potrf(&uplo, &intn, data_a, &intn, &info);

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            return 1;
        }

        if (overwrite_a && clean) {
            zero_other_triangle(uplo, data_a, n);
        }

    } // end of batching loop

    return 1;
}
