/*
 * Templated loops for `vector_norm` and `matrix_norm`
 */
#pragma once

namespace sp_linalg {


// Which matrix norm to compute; python side maps to these values
enum MATRIX_NORM_ORD : Py_ssize_t
{
    FRO = 0,
    ONE = 1,
    INF = 2,
};


template<typename T> int _vector_norm(PyArrayObject *ap_av, PyArrayObject *ap_res, int ord, SliceStatusVec &vec_status) {
    using real_type = real_of_t<T>;
    int info = 0;

    // -------------------------------------------------------------------
    // Input array attributes
    // -------------------------------------------------------------------
    char *A_data = PyArray_BYTES(ap_av);
    int ndim = PyArray_NDIM(ap_av);
    npy_intp *strides = PyArray_STRIDES(ap_av);
    npy_intp *shape = PyArray_DIMS(ap_av);

    real_type *res_data = (real_type *)PyArray_DATA(ap_res);

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 1; i++) { outer_size *= shape[i]; }

    // -------------------------------------------------------------------
    // Batching loop
    // -------------------------------------------------------------------
    T *slice_ptr_a = NULL;
    CBLAS_INT N = shape[ndim-1];
    CBLAS_INT incx = strides[ndim-1] / sizeof(T);
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        // The slices are vectors and `compute_slice_ptr` assumes matrices
        slice_ptr_a = compute_slice_ptr<T>(idx, A_data, ndim + 1, shape, strides);
        res_data[idx] = blas::nrm2(N, slice_ptr_a, incx);
    }

    return info;
}



template <typename T> int _matrix_norm(PyArrayObject *ap_Am, PyArrayObject *ap_res, MATRIX_NORM_ORD ord) {
    using real_type = real_of_t<T>;
    int info = 0;

    // -------------------------------------------------------------------
    // Input array attributes
    // -------------------------------------------------------------------
    char *A_data = PyArray_BYTES(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);

    real_type *res_data = (real_type *)PyArray_DATA(ap_res);

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    char norm = (ord == MATRIX_NORM_ORD::FRO) ? 'F' : (ord == MATRIX_NORM_ORD::ONE) ? '1' : 'I';

    // -------------------------------------------------------------------
    // Workspace allocation
    // -------------------------------------------------------------------
    CBLAS_INT intm = shape[ndim - 2], intn = shape[ndim - 1];
    CBLAS_INT lwork = (ord == MATRIX_NORM_ORD::INF) ? intm : 0;
    CBLAS_INT buff_size = intm * intn;

    real_type *work = (real_type *)PyMem_RawMalloc(lwork * sizeof(real_type));
    if (work == NULL) { return -101; }

    T *buff = (T *)PyMem_RawMalloc(buff_size * sizeof(T));
    if (buff == NULL) { PyMem_FREE(work); return -102; }


    // -------------------------------------------------------------------
    // Main batching loop
    // -------------------------------------------------------------------
    T *slice_ptr_A = NULL;
    real_type *res_ptr = res_data;

    for (npy_intp idx = 0; idx < outer_size; idx++) {
        slice_ptr_A = compute_slice_ptr<T>(idx, A_data, ndim, shape, strides);

        copy_strided_to_f(buff, slice_ptr_A, intm, intn, strides[ndim - 2], strides[ndim - 1], intm);
        *res_ptr = lapack::lange(norm, intm, intn, buff, intm, work);

        // `res_ptr` is C-ordered by construction and just contains scalars.
        res_ptr += 1;
    }

    free(buff);
    free(work);
    return info;
}



} // namespace sp_linalg
