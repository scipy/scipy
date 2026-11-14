/*
 * Templated loops for `vector_norm` and `matrix_norm`
 */
#pragma once

namespace sp_linalg {


template <typename T, bool neg> inline real_of_t<T> _vector_norm_1(T *data, const CBLAS_INT n, const CBLAS_INT incx) {
    real_of_t<T> res = 0.0;

    for (CBLAS_INT i = 0; i < n; i++) {
        T data_ij = data[i * incx];

        if constexpr (neg) {
            res += 1 / std::abs(data_ij);
        } else {
            res += std::abs(data_ij);
        }

    }

    if constexpr (neg) {
        return 1 / res;
    } else {
        return res;
    }
}

template <typename T, bool neg> inline real_of_t<T> _vector_norm_2(T *data, const CBLAS_INT n, const CBLAS_INT incx) {
    real_of_t<T> res = 0.0;

    for (CBLAS_INT i = 0; i < n; i++) {
        T data_ij = data[i * incx];
        real_of_t<T> data_ij_squared;

        if constexpr (is_complex_v<T>) {
            data_ij_squared = std::real(data_ij) * std::real(data_ij) + std::imag(data_ij) * std::imag(data_ij);
        } else {
            data_ij_squared = data_ij * data_ij;
        }

        if constexpr (neg) {
            res += 1 / data_ij_squared;
        } else {
            res += data_ij_squared;
        }

    }

    if constexpr (neg) {
        return 1 / std::sqrt(res);
    } else {
        return std::sqrt(res);
    }
}

/*
 * @brief compute the Frobenius norm of a matrix.
 *
 * @param data  Pointer to the data.
 * @param m     The number of rows in the matrix.
 * @param n     The number of columns in the matrix.
 * @param s2    The strides along the row dimension, in elements.
 * @param s1    The strides along the column dimension, in elements.
 * @param neg   Whether to compute the positive or negative norm.
 */
template <typename T, bool neg> inline real_of_t<T> _matrix_norm_F(T *data, const CBLAS_INT m, const CBLAS_INT n, const CBLAS_INT s2, const CBLAS_INT s1) {
    real_of_t<T> res = 0.0;

    for (CBLAS_INT i = 0; i < n; i++) {
        for (CBLAS_INT j = 0; j < m; j++) {
            T data_ij = data[i * s1 + j * s2];
            real_of_t<T> data_ij_squared;

            if constexpr (is_complex_v<T>) {
                data_ij_squared = std::real(data_ij) * std::real(data_ij) + std::imag(data_ij) * std::imag(data_ij);
            } else {
                data_ij_squared = data_ij * data_ij;
            }

            if constexpr (neg) {
                res += 1 / data_ij_squared;
            } else {
                res += data_ij_squared;
            }
        }
    }

    return std::sqrt(res);
}


/**
 * @brief compute the +-1 norm of a matrix. This can also be used for the infinity norm by swapping row and column dimensions and strides.
 *
 * @param data  Pointer to the data.
 * @param m     The number of rows in the matrix.
 * @param n     The number of columns.
 * @param s2    Strides along the row dimension in elements.
 * @param s1    Strides along the column dimension in elements.
 * @param neg   Whether the compute the +1 or -1 norm.
 */
template <typename T, bool neg> inline real_of_t<T> _matrix_norm_1(T *data, const CBLAS_INT m, const CBLAS_INT n, const CBLAS_INT s2, const CBLAS_INT s1) {
    real_of_t<T> res;
    if constexpr (neg) {
        res = std::numeric_limits<real_of_t<T>>::infinity();
    } else {
        res = 0.0;
    }

    for (CBLAS_INT i = 0; i < n; i++) {
        real_of_t<T> tmp = 0.0;

        for (CBLAS_INT j = 0; j < m; j++) {
            tmp += std::abs(data[i * s1 + j * s2]);
        }

        if constexpr (neg) {
            if (tmp < res) { res = tmp; }
        } else {
            if (tmp > res) { res = tmp; }
        }
    }

    return res;
}


template <typename T> int _vector_norm(PyArrayObject *ap_av, PyArrayObject *ap_res, int ord, SliceStatusVec &vec_status) {
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

        switch (ord) {
            case Vector_norm_ord::V_ONE:
                res_data[idx] = _vector_norm_1<T, false>(slice_ptr_a, intn, incx);
                break;
            case Vector_norm_ord::V_NEG_ONE:
                res_data[idx] = _vector_norm_1<T, true>(slice_ptr_a, intn, incx);
                break;
            case Vector_norm_ord::V_TWO:
                res_data[idx] = _vector_norm_2<T, false>(slice_ptr_a, intn, incx);
                break;
            case Vector_norm_ord::V_NEG_TWO:
                res_data[idx] = _vector_norm_2<T, true>(slice_ptr_a, intn, incx);
                break;
            case Vector_norm_ord::V_INF:
                res_data[idx] = _matrix_norm_1<T, false>(slice_ptr_a, 1, intn, 0, incx);
                break;
            case Vector_norm_ord::V_NEG_INF:
                res_data[idx] = _matrix_norm_1<T, true>(slice_ptr_a, 1, intn, 0, incx);
                break;
            default:
                return -102;
        }

    }

    return info;
}


template <typename T> int _matrix_norm(PyArrayObject *ap_Am, PyArrayObject *ap_res, Matrix_norm_ord ord) {
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
            case Matrix_norm_ord::M_FRO:
                res_ptr[idx] = _matrix_norm_F<T, false>(slice_ptr_A, intm, intn, s2, s1);
                break;
            case Matrix_norm_ord::M_ONE:
                res_ptr[idx] = _matrix_norm_1<T, false>(slice_ptr_A, intm, intn, s2, s1);
                break;
            case Matrix_norm_ord::M_NEG_ONE:
                res_ptr[idx] = _matrix_norm_1<T, true>(slice_ptr_A, intm, intn, s2, s1);
                break;
            case Matrix_norm_ord::M_INF:
                // Effectively transpose to iterate over rows instead of columns
                res_ptr[idx] = _matrix_norm_1<T, false>(slice_ptr_A, intn, intm, s1, s2);
                break;
            case Matrix_norm_ord::M_NEG_INF:
                res_ptr[idx] = _matrix_norm_1<T, true>(slice_ptr_A, intn, intm, s1, s2);
                break;
            default:
                return -102;
        }

    }

    return info;
}



} // namespace sp_linalg
