
namespace sp_linalg {

template<typename T>
int
_cholesky(PyArrayObject *ap_Am, PyArrayObject *ap_Cm, int lower, int overwrite_a, int clean, SliceStatusVec &vec_status) {
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
    npy_intp s0 = strides[ndim-2] / sz;
    npy_intp s1 = strides[ndim-1] / sz;
    int f_contig = (s0 == 1 && s1 == n);

    /*
     * LAPACK's potrf interprets data as F-ordered, but scipy's convention
     * is to return C-ordered output. By flipping uplo, potrf's result lands
     * in the correct C-order triangle directly without any further work:
     *
     * Examples:
     *   - Real C-ordered, lower=True, overwrite_a=False
     *.      1. Lower triangle of A is copied with copy_triangle_to_C
     *.      2. This creates a lower triangular C-ordered array, which is
     *          which is equivalent to an UPPER triangular F-ordered array.
     *       3. potrf only reads from triangle given by uplo and assumes symmetry
     *       4. potrf is called with uplo = "U" overwriting the upper-triangle
     *          of the F-ordered array with U.
     *       5. Scipy sees this as a LOWER triangular C-ordered array and returns it.
     *
     *.  - Complex C-ordered, lower=False, overwrite_a=False
     *.      1. Upper triangle of A is copied with copy_triangle_to_C
     *       2. This creates a LOWER triangular F-ordered array with values
     *          equal to the COMPLEX CONJUGATE conj(A) of the lower triangle of the original
     *          matrix.
     *       3. potrf only reads from triangle given by uplo and assumes Hermitian
     *       4. potrf is called with uplo = "L" overwriting the lower-triangle of
     *          F-ordered array Am* with its lower Cholesky decomposition conj(L).
     *       5. Scipy sees this as an UPPER triangular C-ordered array U = conj(L).T = L^H,
     *          which is exactly the relation between the lower and upper Cholesky
     *          decompositions and this can be returned directly to the user.
     *
     *  - Real/Complex, F-ordered, lower=True, overwrite_a=True
     *      1. potrf is called with uplo = "L" overwriting the lower-triangle of
     *         the F-ordered array with L* and this is returned directly to the user.
     *     This is the only case where scipy returns F-ordered output.
     */
    char uplo_f = (overwrite_a && f_contig) ? (lower ? 'L' : 'U')
                                             : (lower ? 'U' : 'L');

    T *data_a = NULL;

    // Main loop over batch
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        // Ensure data is in the correct triangle of the F-ordered input data_a
        if (!overwrite_a) {
            data_a = &ret_data[idx * n * n];
            T *slice_ptr_A = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_triangle_to_C(data_a, slice_ptr_A, n, n, s0, s1, lower ? 'L' : 'U');
        } else {
            data_a = Am_data;
        }
        // NB. `overwrite_a` is only enabled when the input is 2D contiguous so
        // data is already correctly ordered in input array. If generalized to
        // `ndim > 2` this should be updated accordingly.

        init_status(slice_status, idx, slice_structure);
        call_potrf(&uplo_f, &intn, data_a, &intn, &info);

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            return 1;
        }

        if (overwrite_a && clean) {
            zero_other_triangle(uplo_f, data_a, n);
        }

    } // end of batching loop

    return 1;
}

} // namespace sp_linalg
