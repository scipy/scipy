#include "npy_cblas.h"
#include "_npymath.hh"
#include "_common_array_utils.hh"


/*
 * SVD size helper: if A.shape == (m, n),
 * U is either (m, m) or (m, k) and Vh is either (n, n) or (k, n)
 */ 
int
u_vh_shapes(npy_intp m, npy_intp n, char jobz,
            npy_intp *u_shape0, npy_intp *u_shape1,
            npy_intp *vh_shape0, npy_intp *vh_shape1
){
    npy_intp k = m < n ? m : n; 

    switch(jobz) {
        case('N') :
            *u_shape0 = 0;
            *u_shape1 = 0;
            *vh_shape0 = 0;
            *vh_shape1 = 0;
            break;
        case('A'):
            *u_shape0 = m;
            *u_shape1 = m;
            *vh_shape0 = n;
            *vh_shape1 = n;
            break;
        case('S'):
            *u_shape0 = m;
            *u_shape1 = k;
            *vh_shape0 = k;
            *vh_shape1 = n;
            break;
        default:
            return 0;  // error
    }
    return 1;
}


template<typename T>
int
_svd_gesdd(PyArrayObject* ap_Am, PyArrayObject *ap_U, PyArrayObject *ap_S, PyArrayObject *ap_Vh, char jobz, SliceStatusVec& vec_status)
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
    if (ndim > 2) {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    // output array attributes
    npy_intp min_mn = m < n ? m : n;
    npy_intp max_mn = m < n ? n : m;

    // U, Vh shapes ( =0 if compute_uv=False)
    npy_intp u_shape0, u_shape1, vh_shape0, vh_shape1;
    if(!u_vh_shapes(m, n, jobz, &u_shape0, &u_shape1, &vh_shape0, &vh_shape1)) {
        return -101; // unknown value for jobz, bail out
    }

    CBLAS_INT ldu = (CBLAS_INT)(u_shape0 == 0 ? 1 : u_shape0);
    CBLAS_INT ldvh = (CBLAS_INT)(vh_shape0 == 0 ? 1 : vh_shape0);

    T *ptr_U = u_shape0 > 0 ? (T *)PyArray_DATA(ap_U) : NULL;
    T *ptr_Vh = vh_shape0 > 0 ? (T *)PyArray_DATA(ap_Vh) : NULL;
    real_type *ptr_S = (real_type *)PyArray_DATA(ap_S);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, intm = (CBLAS_INT)m, lwork = -1, info;
    T tmp = numeric_limits<T>::zero;

    // query LWORK
    gesdd(&jobz, &intm, &intn, NULL, &intm, NULL, NULL, &ldu, NULL, &ldvh, &tmp, &lwork, NULL, NULL, &info);
    if (info != 0) { info = -100; return (int)info; }

    lwork = (CBLAS_INT)(real_part(tmp));
    if(lwork == 0) { lwork = 1; }

    // allocate
    npy_intp bufsize = m*n + lwork;
    if (jobz != 'N') {
        bufsize += u_shape0 * u_shape1 + vh_shape0 * vh_shape1;    // U and Vh, if referenced
    }

    T *buf = (T *)malloc(bufsize*sizeof(T));
    if (buf == NULL) { info = -101; return (int)info; }

    // partition the workspace
    T *data = &buf[0];
    T *work = &buf[m*n];

    T *buf_U = NULL;
    T *buf_Vh = NULL;
    if (jobz != 'N') {
        buf_U = &buf[m*n + lwork];
        buf_Vh = &buf[m*n + lwork + u_shape0*u_shape1];
    }

    CBLAS_INT *iwork = NULL;
    real_type *rwork = NULL;
    // iwork
    iwork = (CBLAS_INT *)malloc(8*min_mn*sizeof(CBLAS_INT));
    if (iwork == NULL) {
        free(buf);
        info = -102;
        return (int)info;
    }

    // rwork
    if (type_traits<T>::is_complex) {
        // assume LAPACK > 3.6 (cf LAPACK docs on netlib.org)
        npy_intp lrwork = std::max(
            5*min_mn*min_mn + 5*min_mn,
            2*max_mn*min_mn + 2*min_mn*min_mn + min_mn
        );
        rwork = (real_type *)malloc(lrwork * sizeof(real_type));
        if (rwork == NULL) {
            free(buf);
            free(iwork);
            info = -103;
            return (int)info;
        }
    }

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        // copy the slice to `data` in F order
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
        copy_slice_F(data, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);

        // SVD the slice
        gesdd(&jobz, &intm, &intn, data, &intm, ptr_S, buf_U, &ldu, buf_Vh, &ldvh, work, &lwork, rwork, iwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy-and-tranpose U and Vh slices from temp buffers to the output;
        // S slice is filled in in-place already; 
        copy_slice_F_to_C(ptr_U, buf_U, u_shape0, u_shape1);
        copy_slice_F_to_C(ptr_Vh, buf_Vh, vh_shape0, vh_shape1);

        // advance the output pointers: U, S, Vh are C-contiguous by construction
        ptr_U += u_shape0 * u_shape1;
        ptr_Vh += vh_shape0 * vh_shape1;
        ptr_S += min_mn;
    }

 done:
    free(buf);
    free(iwork);
    free(rwork);
    return 0;
}


template<typename T>
int
_svd_gesvd(PyArrayObject* ap_Am, PyArrayObject *ap_U, PyArrayObject *ap_S, PyArrayObject *ap_Vh, char jobz, SliceStatusVec& vec_status)
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
    if (ndim > 2) {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    // output array attributes
    npy_intp min_mn = m < n ? m : n;

    // U, Vh shapes ( =0 if compute_uv=False)
    npy_intp u_shape0, u_shape1, vh_shape0, vh_shape1;
    if(!u_vh_shapes(m, n, jobz, &u_shape0, &u_shape1, &vh_shape0, &vh_shape1)) {
        return -101; // unknown value for jobz, bail out
    }

    CBLAS_INT ldu = (CBLAS_INT)(u_shape0 == 0 ? 1 : u_shape0);
    CBLAS_INT ldvh = (CBLAS_INT)(vh_shape0 == 0 ? 1 : vh_shape0);

    T *ptr_U = u_shape0 > 0 ? (T *)PyArray_DATA(ap_U) : NULL;
    T *ptr_Vh = vh_shape0 > 0 ? (T *)PyArray_DATA(ap_Vh) : NULL;
    real_type *ptr_S = (real_type *)PyArray_DATA(ap_S);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn = (CBLAS_INT)n, intm = (CBLAS_INT)m, lwork = -1, info;
    T tmp = numeric_limits<T>::zero;

    // query LWORK
    gesvd(&jobz, &jobz, &intm, &intn, NULL, &intm, NULL, NULL, &ldu, NULL, &ldvh, &tmp, &lwork, NULL, &info);
    if (info != 0) { info = -100; return (int)info; }

    lwork = (CBLAS_INT)(real_part(tmp));
    if(lwork == 0) { lwork = 1; }

    // allocate
    npy_intp bufsize = m*n + lwork;
    if (jobz != 'N') {
        bufsize += u_shape0 * u_shape1 + vh_shape0 * vh_shape1;    // U and Vh, if referenced
    }

    T *buf = (T *)malloc(bufsize*sizeof(T));
    if (buf == NULL) { info = -101; return (int)info; }

    // partition the workspace
    T *data = &buf[0];
    T *work = &buf[m*n];

    T *buf_U = NULL;
    T *buf_Vh = NULL;
    if (jobz != 'N') {
        buf_U = &buf[m*n + lwork];
        buf_Vh = &buf[m*n + lwork + u_shape0*u_shape1];
    }

    CBLAS_INT *iwork = NULL;
    real_type *rwork = NULL;
    rwork = (real_type *)malloc(5*min_mn*sizeof(real_type));
    if (rwork == NULL) {
        free(buf);
        info = -103;
        return (int)info;
    }

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        // copy the slice to `data` in F order
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
        copy_slice_F(data, slice_ptr, m, n, strides[ndim-2], strides[ndim-1]);

        // SVD the slice
        gesvd(&jobz, &jobz, &intm, &intn, data, &intm, ptr_S, buf_U, &ldu, buf_Vh, &ldvh, work, &lwork, rwork, &info);

        if(info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            goto done;
        }

        // copy-and-tranpose U and Vh slices from temp buffers to the output;
        // S slice is filled in in-place already; 
        copy_slice_F_to_C(ptr_U, buf_U, u_shape0, u_shape1);
        copy_slice_F_to_C(ptr_Vh, buf_Vh, vh_shape0, vh_shape1);

        // advance the output pointers: U, S, Vh are C-contiguous by construction
        ptr_U += u_shape0 * u_shape1;
        ptr_Vh += vh_shape0 * vh_shape1;
        ptr_S += min_mn;
    }

 done:
    free(buf);
    free(iwork);
    free(rwork);
    return 0;
}


template<typename T>
int
_svd(PyArrayObject* ap_Am, PyArrayObject *ap_U, PyArrayObject *ap_S, PyArrayObject *ap_Vh, char jobz, const char * lapack_driver, SliceStatusVec& vec_status)
{
    int info;
    if (strcmp(lapack_driver, "gesdd") == 0) {
        info = _svd_gesdd<T>(ap_Am, ap_U, ap_S, ap_Vh, jobz, vec_status);
    }
    else if (strcmp(lapack_driver, "gesvd") == 0) {
        info = _svd_gesvd<T>(ap_Am, ap_U, ap_S, ap_Vh, jobz, vec_status);
    }
    else {
        // should have been validated at call site, really
        info = -110;
    }
    return info;
}
