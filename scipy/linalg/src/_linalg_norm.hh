/*
 * Templated loops for `vector_norm` and `matrix_norm`
 */
#pragma once

namespace sp_linalg {



/*
 * `s2` denote strides along the rows, `s1` denote strides along the columns.
 */
template<typename T> inline real_of_t<T> _matrix_norm_F(T *data, const CBLAS_INT m, const CBLAS_INT n, const CBLAS_INT s2, const CBLAS_INT s1) {
    real_of_t<T> res = 0.0;

    for (CBLAS_INT i = 0; i < n; i++) {
        for (CBLAS_INT j = 0; j < m; j++) {
            T data_ij = data[i * s1 + j * s2];
            if constexpr (is_complex_v<T>) {
                res += std::real(data_ij) * std::real(data_ij) + std::imag(data_ij) * std::imag(data_ij);
            } else {
                res += data_ij * data_ij;
            }
        }
    }

    return std::sqrt(res);
}


template<typename T> inline real_of_t<T> _matrix_norm_1(T *data, CBLAS_INT m, CBLAS_INT n, const CBLAS_INT s2, const CBLAS_INT s1) {
    real_of_t<T> res = 0.0;

    for (CBLAS_INT i = 0; i < n; i++) {
        real_of_t<T> tmp = 0.0;

        for (CBLAS_INT j = 0; j < m; j++) {
            tmp += std::abs(data[i * s1 + j * s2]);
        }

        if (tmp > res) { res = tmp; }
    }

    return res;
}


template<typename T> int _vector_norm(PyArrayObject *ap_av, PyArrayObject *ap_res, int ord, SliceStatusVec &vec_status) {
    using real_type = real_of_t<T>;
    int info = 0;

    // -------------------------------------------------------------------
    // Input array attributes
    // -------------------------------------------------------------------
    T *A_data = (T *)PyArray_DATA(ap_av);
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
    CBLAS_INT intn = shape[ndim-1];
    CBLAS_INT incx = strides[ndim-1] / sizeof(T);
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        // The slices are vectors and `compute_slice_ptr` assumes matrices
        slice_ptr_a = compute_slice_ptr(idx, A_data, ndim + 1, shape, strides);
        res_data[idx] = _matrix_norm_F(slice_ptr_a, intn, 1, incx, 0);
    }

    return info;
}



template <typename T> int _matrix_norm(PyArrayObject *ap_Am, PyArrayObject *ap_res, MATRIX_NORM_ORD ord) {
    using real_type = real_of_t<T>;
    int info = 0;

    // -------------------------------------------------------------------
    // Input array attributes
    // -------------------------------------------------------------------
    T *A_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);

    real_type *res_data = (real_type *)PyArray_DATA(ap_res);

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    // -------------------------------------------------------------------
    // Main batching loop
    // -------------------------------------------------------------------
    T *slice_ptr_A = NULL;
    real_type *res_ptr = res_data;
    CBLAS_INT intm = shape[ndim - 2], intn = shape[ndim - 1];

    CBLAS_INT s2 = strides[ndim - 2] / sizeof(T);
    CBLAS_INT s1 = strides[ndim - 1] / sizeof(T);

    for (npy_intp idx = 0; idx < outer_size; idx++) {
        slice_ptr_A = compute_slice_ptr(idx, A_data, ndim, shape, strides);

        switch (ord) {
            case MATRIX_NORM_ORD::FRO:
                *res_ptr = _matrix_norm_F(slice_ptr_A, intm, intn, s2, s1);
                break;
            case MATRIX_NORM_ORD::ONE:
                *res_ptr = _matrix_norm_1(slice_ptr_A, intm, intn, s2, s1);
                break;
            case MATRIX_NORM_ORD::INF:
                // Effectively transpose to iterate over rows instead of columns
                *res_ptr = _matrix_norm_1(slice_ptr_A, intn, intm, s1, s2);
                break;
            default:
                return -102;
        }

        // `res_ptr` is C-ordered by construction and just contains scalars.
        res_ptr += 1;
    }

    return info;
}



} // namespace sp_linalg
