/*
 * Templated loops for linalg.eig
 */


/*
 * Copy eigenvectors from F-ordered real `v` (ldv, n)-shaped array in the ?GEEV convention
 * into an (n, n)-shaped C-ordered complex array `dst`.
 */
template<typename real_type, typename cmplx_type>
void
transform_eigvecs(cmplx_type dst, real_type *v, CBLAS_INT ldv, CBLAS_INT n, real_type *wi) {
    for (CBLAS_INT j=0; j < n; ) {
        if(wi[j] == 0.) {
            // If the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column of VL.
            for (CBLAS_INT i=0; i<n; i++) {
                dst[i*n + j] = cpack(v[i + j*ldv], real_type(0.));
            }
            j += 1;
        }
        else {
            // If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
            // then u(j) = VL(:,j) + i*VL(:,j+1) and u(j+1) = VL(:,j) - i*VL(:,j+1).
            for (CBLAS_INT i=0; i<n; i++) {
                real_type re = v[i + j*ldv], im = v[i + (j+1)*ldv];
                dst[i*n + j] = cpack(re, im);     // VL(i, j)
                dst[i*n + j + 1] = cpack(re, -im);  // VL(i, j+1)
            }
            j += 2;
        }
    }
}


template<typename T>
int
_eig(PyArrayObject* ap_Am, PyArrayObject *ap_w, PyArrayObject *ap_vl, PyArrayObject *ap_vr, SliceStatusVec& vec_status)
{
    using real_type = typename type_traits<T>::real_type; // float if T==npy_cfloat etc
    using npy_complex_type = typename type_traits<T>::npy_complex_type;
    SliceStatus slice_status;

    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    T* Am_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];
    npy_intp* strides = PyArray_STRIDES(ap_Am);
    
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2) {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    // Output array pointers
    npy_complex_type *ptr_W = (npy_complex_type *)PyArray_DATA(ap_w);

    int compute_vl = (ap_vl != NULL);
    int compute_vr = (ap_vr != NULL);

    npy_complex_type *ptr_vl = compute_vl ? (npy_complex_type *)PyArray_DATA(ap_vl) : NULL;
    npy_complex_type *ptr_vr = compute_vr ? (npy_complex_type *)PyArray_DATA(ap_vr) : NULL;

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, lwork = -1, info;
    T tmp = numeric_limits<T>::zero;

    char jobvl = compute_vl ? 'V': 'N', jobvr = compute_vr ? 'V' : 'N';
    CBLAS_INT lda = n;
    CBLAS_INT ldvl = n;
    CBLAS_INT ldvr = n;

    // c- and z variants: lwork query segfaults with rwork=NULL, allocate it straight away
    real_type *rwork = NULL;
    if (type_traits<T>::is_complex) {
        rwork = (real_type *)malloc(2*n*sizeof(real_type));
        if (rwork == NULL) {
            return -100;
        }
    }

    // query LWORK
    geev(&jobvr, &jobvl, &intn, NULL, &lda, NULL, NULL, NULL, &ldvl, NULL, &ldvr, &tmp, &lwork, rwork, &info);
    if (info != 0) { free(rwork);  return -101; }

    lwork = _calc_lwork(tmp);
    if (lwork < 0) { free(rwork); return -102; }

    // allocate
    npy_intp bufsize = n*n + lwork + n;
    npy_intp wi_size = type_traits<T>::is_complex ? 0 : n ;
    bufsize += wi_size;

    npy_intp vl_size = compute_vl ? ldvl*n : 0;
    npy_intp vr_size = compute_vr ? ldvr*n : 0;
    bufsize += vl_size + vr_size;

    T *buf = (T *)malloc(bufsize*sizeof(T));
    if (buf == NULL) { free(rwork); return -103; }

    // partition the workspace
    T *data = &buf[0];
    T *work = &buf[n*n];
    T *wr = &buf[n*n + lwork];

    T *wi = NULL;
    if(wi_size > 0) { wi = &buf[n*n + lwork + n]; }

    T *buf_vl = NULL;
    if(vl_size > 0) { buf_vl = &buf[n*n + lwork + n + wi_size]; }

    T *buf_vr = NULL;
    if(vr_size > 0) { buf_vr = &buf[n*n + lwork + n + wi_size + vl_size]; }

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        // copy the slice to `data` in F order
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
        copy_slice_F(data, slice_ptr, n, n, strides[ndim-2], strides[ndim-1]);

        // compute eigenvalues for the slice
        geev(&jobvl, &jobvr, &intn, data, &lda, wr, wi, buf_vl, &ldvl, buf_vr, &ldvr, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy-and-tranpose W, VR and VL slices from temp buffers to the output;
        if constexpr (type_traits<T>::is_complex) {
            memcpy(ptr_W, wr, n*sizeof(T));
            ptr_W += n;

            if (compute_vl) {
                copy_slice_F_to_C(ptr_vl, buf_vl, n, n, ldvl);
                ptr_vl += n*n;
            }
            if (compute_vr) {
                copy_slice_F_to_C(ptr_vr, buf_vr, n, n, ldvr);
                ptr_vr += n*n;
            }
        }
        else {
            // convert wr,wi into w
            for(npy_intp i=0; i<n; i++) {
                ptr_W[i] = cpack(wr[i], wi[i]);
            }
            ptr_W += n;

            if (compute_vl) {
                transform_eigvecs(ptr_vl, buf_vl, ldvl, n, wi);
                ptr_vl += n*n;
            }
            if (compute_vr) {
                transform_eigvecs(ptr_vr, buf_vr, ldvr, n, wi);
                ptr_vr += n*n;
            }
        }
    }

 done:
    free(buf);
    free(rwork);

    return 0;
}
