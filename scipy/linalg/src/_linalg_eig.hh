#pragma once
/*
 * Templated loops for linalg.eig
 */


namespace sp_linalg {

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
                dst[i*n + j] = detail::cpack(v[i + j*ldv], real_type(0.));
            }
            j += 1;
        }
        else {
            // If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
            // then u(j) = VL(:,j) + i*VL(:,j+1) and u(j+1) = VL(:,j) - i*VL(:,j+1).
            for (CBLAS_INT i=0; i<n; i++) {
                real_type re = v[i + j*ldv], im = v[i + (j+1)*ldv];
                dst[i*n + j] = detail::cpack(re, im);     // VL(i, j)
                dst[i*n + j + 1] = detail::cpack(re, -im);  // VL(i, j+1)
            }
            j += 2;
        }
    }
}


template<typename T>
int
_reg_eig(PyArrayObject* ap_Am, PyArrayObject *ap_w, PyArrayObject *ap_vl, PyArrayObject *ap_vr, int overwrite_a, SliceStatusVec& vec_status)
{
    using real_type = typename detail::type_traits<T>::real_type; // float if T==npy_cfloat etc
    using npy_complex_type = typename detail::type_traits<T>::npy_complex_type;
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
    T tmp = detail::numeric_limits<T>::zero;

    char jobvl = compute_vl ? 'V': 'N', jobvr = compute_vr ? 'V' : 'N';
    CBLAS_INT lda = n;
    CBLAS_INT ldvl = n;
    CBLAS_INT ldvr = n;

    // c- and z variants: lwork query segfaults with rwork=NULL, allocate it straight away
    real_type *rwork = NULL;
    if constexpr (detail::type_traits<T>::is_complex) {
        rwork = (real_type *)malloc(2*n*sizeof(real_type));
        if (rwork == NULL) {
            return -100;
        }
    }

    // query LWORK
    call_geev(&jobvl, &jobvr, &intn, NULL, &lda, NULL, NULL, NULL, &ldvl, NULL, &ldvr, &tmp, &lwork, rwork, &info);
    if (info != 0) { free(rwork);  return -101; }

    lwork = _calc_lwork(tmp);
    if (lwork < 0) { free(rwork); return -102; }

    /*
     * Allocate memory and chop the buffer into parts
     *
     *     lwork        n      data_size   wi_size
     * |-----------|---------|-----------|----------|---------|-------|
     * ^           ^         ^           ^          ^         ^
     * work        wr        data        wi         buf_vl    buf_vr
     *
     * Here
     *   - data is n*n if overwrite_a is False (and =0 otherwise)
     *   - wr & wi are eigenvalues;
     *     wr is always length n; wi is only used for real inputs (s- and d-geev
     *     have both in `wr` and `wi` arguments, while c- and z-geev only have `wr`)
     *   - `buf_vl` and `buf_vr` hold eigenvectors if requested via `compute_{vl,vr}`.
     *
     * NB: we do not implement jobz='O' yet, so we never reuse A for U or Vh.
     */
    npy_intp data_size = overwrite_a ? 0 : n*n;
    npy_intp wi_size = detail::type_traits<T>::is_complex ? 0 : n;
    npy_intp bufsize = data_size + wi_size + lwork + n;

    npy_intp vl_size = compute_vl ? ldvl*n : 0;
    npy_intp vr_size = compute_vr ? ldvr*n : 0;
    bufsize += vl_size + vr_size;

    T *buf = (T *)malloc(bufsize*sizeof(T));
    if (buf == NULL) { free(rwork); return -103; }

    // partition the workspace
    T *work = &buf[0];
    T *wr = &buf[lwork];

    T *data = NULL;
    if (overwrite_a) {
        // work in-place (2D only, ensured at the python call site)
        data = (T *)Am_data;
    } else {
        data = &buf[lwork + n];
    }

    T *wi = NULL;
    if(wi_size > 0) { wi = &buf[lwork + n + data_size]; }

    T *buf_vl = NULL;
    if(vl_size > 0) { buf_vl = &buf[lwork + n + data_size + wi_size]; }

    T *buf_vr = NULL;
    if(vr_size > 0) { buf_vr = &buf[lwork + n + data_size + wi_size + vl_size]; }

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            // copy the slice to `data` in F order
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_slice_F(data, slice_ptr, n, n, strides[ndim-2], strides[ndim-1]);
        }
        // NB: overwrite_a is for 2D F-ordered input only; if it is ever generalized to ndim>2,
        // will need to adjust the data pointer here, too.

        // compute eigenvalues for the slice
        call_geev(&jobvl, &jobvr, &intn, data, &lda, wr, wi, buf_vl, &ldvl, buf_vr, &ldvr, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy-and-tranpose W, VR and VL slices from temp buffers to the output;
        if constexpr (detail::type_traits<T>::is_complex) {
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
                ptr_W[i] = detail::cpack(wr[i], wi[i]);
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


template<typename T>
int
_gen_eig(PyArrayObject* ap_Am, PyArrayObject *ap_Bm, PyArrayObject *ap_w, PyArrayObject *ap_beta, PyArrayObject *ap_vl, PyArrayObject *ap_vr, int overwrite_a, int overwrite_b, SliceStatusVec& vec_status)
{
    using real_type = typename detail::type_traits<T>::real_type; // float if T==npy_cfloat etc
    using npy_complex_type = typename detail::type_traits<T>::npy_complex_type;
    SliceStatus slice_status;

    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    T* Am_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];
    npy_intp *strides_A = PyArray_STRIDES(ap_Am);

    T *Bm_data = (T *)PyArray_DATA(ap_Bm);
    // shape(B) == shape(A)
    npy_intp *strides_B = PyArray_STRIDES(ap_Bm);

    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2) {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    // Output array pointers
    npy_complex_type *ptr_W = (npy_complex_type *)PyArray_DATA(ap_w);
    T *ptr_beta = (T *)PyArray_DATA(ap_beta);

    int compute_vl = (ap_vl != NULL);
    int compute_vr = (ap_vr != NULL);

    npy_complex_type *ptr_vl = compute_vl ? (npy_complex_type *)PyArray_DATA(ap_vl) : NULL;
    npy_complex_type *ptr_vr = compute_vr ? (npy_complex_type *)PyArray_DATA(ap_vr) : NULL;

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, lwork = -1, info;
    T tmp = detail::numeric_limits<T>::zero;

    char jobvl = compute_vl ? 'V': 'N', jobvr = compute_vr ? 'V' : 'N';
    CBLAS_INT lda = n;
    CBLAS_INT ldb = n;
    CBLAS_INT ldvl = n;
    CBLAS_INT ldvr = n;

    // similar to geev, allocate rwork right away (not sure if ?ggev segfaults otherwise, too)
    real_type *rwork = NULL;
    if constexpr (detail::type_traits<T>::is_complex) {
        rwork = (real_type *)malloc(8*n*sizeof(real_type));
        if (rwork == NULL) {
            return -100;
        }
    }

    // query LWORK
    call_ggev(&jobvl, &jobvr, &intn, NULL, &lda, NULL, &ldb, NULL, NULL, NULL, NULL, &ldvl, NULL, &ldvr, &tmp, &lwork, rwork, &info);
    if (info != 0) { free(rwork);  return -101; }

    lwork = _calc_lwork(tmp);
    if (lwork < 0) { free(rwork); return -102; }

    /*
     * Allocate memory and chop the buffer into parts
     *
     *     lwork        n        n        n/0      n*n/0    n*n/0    n*n/0   n*n/0
     * |-----------|---------|-------|----------|---------|-------|--------|--------|
     * ^           ^         ^       ^          ^         ^       ^        ^
     * work        alphar    beta    alphai     data_A    data_B  buf_vl   buf_vr
     *
     * Here
     *   - A_size & B_size are n*n if overwrite_{a,b} is False (and =0 otherwise)
     *   - alphar & beta are eigenvalues, size `n` (always allocated)
     *   - alphai is the imaginary part of eigenvalues, size `n` if real, 0 otherwise
     *     (s- and d-ggev have both in `alphar` and `alphai` arguments, while c- and
     *     z-geev only have `alphar`)
     *   - `buf_vl` and `buf_vr` hold eigenvectors if requested via `compute_{vl,vr}`.
     *
     * NB: we do not implement jobz='O' yet, so we never reuse A for U or Vh.
     */
    npy_intp alphai_size = detail::type_traits<T>::is_complex ? 0 : n ;
    npy_intp A_size = overwrite_a ? 0 : n*n;
    npy_intp B_size = overwrite_b ? 0 : n*n;

    npy_intp vl_size = compute_vl ? ldvl*n : 0;
    npy_intp vr_size = compute_vr ? ldvr*n : 0;

    npy_intp bufsize = lwork + n + n + alphai_size + A_size + B_size  + vl_size + vr_size;

    T *buf = (T *)malloc(bufsize*sizeof(T));
    if (buf == NULL) { free(rwork); return -103; }

    T *work = &buf[0];
    T *alphar = &buf[lwork];
    T *beta = &buf[lwork + n];

    T *alphai = NULL;
    if(alphai_size > 0) { alphai = &buf[lwork + 2*n]; }

    T *data_A = NULL;
    if (overwrite_a) {
        // work in-place (2D only, ensured at the python call site)
        data_A = (T *)Am_data;
    } else {
        data_A = &buf[lwork + 2*n + alphai_size];
    }

    T *data_B = NULL;
    if (overwrite_b) {
        // work in-place (2D only, ensured at the python call site)
        data_B = (T *)Bm_data;
    } else {
        data_B = &buf[lwork + 2*n + alphai_size + A_size];
    }

    T *buf_vl = NULL;
    if(vl_size > 0) { buf_vl = &buf[lwork + 2*n + alphai_size + A_size + B_size]; }

    T *buf_vr = NULL;
    if(vr_size > 0) { buf_vr = &buf[lwork + 2*n + alphai_size + vl_size + A_size + B_size]; }


    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            // copy the slice to `data` in F order
            T *slice_ptr_A = compute_slice_ptr(idx, Am_data, ndim, shape, strides_A);
            copy_slice_F(data_A, slice_ptr_A, n, n, strides_A[ndim-2], strides_A[ndim-1]);
        }
        // NB: overwrite_a is for 2D F-ordered input only; if it is ever generalized to ndim>2,
        // will need to adjust the data pointer here, too.

        if (!overwrite_b) {
            T *slice_ptr_B = compute_slice_ptr(idx, Bm_data, ndim, shape, strides_B);
            copy_slice_F(data_B, slice_ptr_B, n, n, strides_B[ndim-2], strides_B[ndim-1]);
        }
        // same deal as with overwrite_a


        // compute eigenvalues for the slice
        call_ggev(&jobvl, &jobvr, &intn, data_A, &lda, data_B, &ldb, alphar, alphai, beta, buf_vl, &ldvl, buf_vr, &ldvr, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy-and-tranpose W, VR and VL slices from temp buffers to the output;
        if constexpr (detail::type_traits<T>::is_complex) {
            // alphar and beta are complex and compatible with the W array 
            memcpy(ptr_W, alphar, n*sizeof(T));
            ptr_W += n;

            memcpy(ptr_beta, beta, n*sizeof(T));
            ptr_beta += n;

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
            // convert alphar,alphai,beta into w
            for(npy_intp i=0; i<n; i++) {
                ptr_W[i] = detail::cpack(alphar[i], alphai[i]);
            }
            ptr_W += n;

            memcpy(ptr_beta, beta, n*sizeof(T));
            ptr_beta += n;

            if (compute_vl) {
                transform_eigvecs(ptr_vl, buf_vl, ldvl, n, alphai);
                ptr_vl += n*n;
            }
            if (compute_vr) {
                transform_eigvecs(ptr_vr, buf_vr, ldvr, n, alphai);
                ptr_vr += n*n;
            }
        }
    }


 done:
    free(buf);
    free(rwork);

    return 0;
}

template<typename T>
int
_eig(PyArrayObject* ap_Am, PyArrayObject *ap_Bm,
     PyArrayObject *ap_w, PyArrayObject *ap_beta,
     PyArrayObject *ap_vl, PyArrayObject *ap_vr,
     int overwrite_a, int overwrite_b,
     SliceStatusVec& vec_status
) {
    int info;
    if (ap_Bm == NULL) {
        // sanity check: B and beta are either both NULL or both not NULL (for a generalized eig problem) 
        if (ap_beta != NULL) { return -222; }

        info = _reg_eig<T>(ap_Am, ap_w, ap_vl, ap_vr, overwrite_a, vec_status);
    }
    else {
        if (ap_beta == NULL) {return -223; }

        info = _gen_eig<T>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
    }
    return info;
}

} // namespace sp_linalg
