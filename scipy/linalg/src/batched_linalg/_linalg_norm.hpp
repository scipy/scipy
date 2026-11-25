/*
 * Templated loops for `vector_norm` and `matrix_norm`
 */
#pragma once

namespace sp_linalg {


// Which matrix norm to compute; python side maps to these values
enum Matrix_norm_ord : Py_ssize_t
{
    m_fro       = 0,
    m_nuc       = 1,
    m_one       = 10,
    m_neg_one   = 11,
    m_two       = 20,
    m_neg_two   = 21,
    m_inf       = 30,
    m_neg_inf   = 31,
};


/*
 * @brief compute the +- 1 norm of a vector `data`
 *
 * @param data  The array containing the data.
 * @param n     The number of elements in the vector.
 * @param incx  The increment between subsequent elements in bytes.
 * @param neg   Whether to compute the + or - 1 norm.
 */
template <typename T, bool neg> inline real_of_t<T> _vector_norm_1(const T *data, const CBLAS_INT n, const CBLAS_INT incx) {
    const char *base = reinterpret_cast<const char *>(data);

    real_of_t<T> res = 0.0;
    for (CBLAS_INT i = 0; i < n; i++) {
        const T data_ij = *reinterpret_cast<const T *>(base + i * incx);

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


/*
 * @brief compute the +- 2 norm of a vector `data`.
 *
 * @param data  The array containing the data.
 * @param n     The number of elements in the vector.
 * @param incx  The increment between subsequent elements in bytes.
 * @param neg   Whether to compute the + or - 2 norm.
 */
template <typename T, bool neg> inline real_of_t<T> _vector_norm_2(const T *data, const CBLAS_INT n, const CBLAS_INT incx) {
    const char *base = reinterpret_cast<const char *>(data);

    real_of_t<T> res = 0.0;
    for (CBLAS_INT i = 0; i < n; i++) {
        const T data_ij = *reinterpret_cast<const T *>(base + i * incx);
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
 * @brief compute the general p norm of vector
 *
 * @param data  Pointer to the data.
 * @param p     The exponent used in the computation of the p norm.
 * @param n     The number of elements in the vector.
 * @param incx  The increment between subsequent elements in the vector, in bytes.
 *
 * @return the p norm
 */
template <typename T> inline real_of_t<T> _vector_norm_general(const T *data, const f64 p, const CBLAS_INT n, const CBLAS_INT incx) {
    const char *base = reinterpret_cast<const char *>(data);

    real_of_t<T> res = 0.0;

    for (CBLAS_INT i = 0; i < n; i++) {
        const T data_ij = *reinterpret_cast<const T *>(base + i * incx);
        res += std::pow(std::abs(data_ij), p);
    }

    return std::pow(res, 1 / p);
}


/*
 * @brief compute the Frobenius norm of a matrix.
 *
 * @param data  Pointer to the data.
 * @param m     The number of rows in the matrix.
 * @param n     The number of columns in the matrix.
 * @param s0    The strides along the row dimension, in bytes.
 * @param s1    The strides along the column dimension, in bytes.
 * @param neg   Whether to compute the positive or negative norm.
 */
template <typename T, bool neg> inline real_of_t<T> _matrix_norm_F(const T *data, const CBLAS_INT m, const CBLAS_INT n, const CBLAS_INT s0, const CBLAS_INT s1) {
    const char *base = reinterpret_cast<const char *>(data);

    real_of_t<T> res = 0.0;
    for (CBLAS_INT i = 0; i < n; i++) {
        for (CBLAS_INT j = 0; j < m; j++) {
            const T data_ij = *reinterpret_cast<const T *>(base + j * s0 + i * s1);
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
 * @param s0    Strides along the row dimension in bytes.
 * @param s1    Strides along the column dimension in bytes.
 * @param neg   Whether the compute the +1 or -1 norm.
 */
template <typename T, bool neg> inline real_of_t<T> _matrix_norm_1(const T *data, const CBLAS_INT m, const CBLAS_INT n, const CBLAS_INT s0, const CBLAS_INT s1) {
    const char *base = reinterpret_cast<const char *>(data);

    real_of_t<T> res;
    if constexpr (neg) {
        res = std::numeric_limits<real_of_t<T>>::infinity();
    } else {
        res = 0.0;
    }

    for (CBLAS_INT i = 0; i < n; i++) {
        real_of_t<T> tmp = 0.0;

        for (CBLAS_INT j = 0; j < m; j++) {
            tmp += std::abs(*reinterpret_cast<const T *>(base + j * s0 + i * s1));
        }

        if constexpr (neg) {
            if (tmp < res) { res = tmp; }
        } else {
            if (tmp > res) { res = tmp; }
        }
    }

    return res;
}


template <typename T> BatchStatus _vector_norm(PyArrayObject *ap_av, PyArrayObject *ap_res, f64 ord) {
    using real_type = real_of_t<T>;

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
    CBLAS_INT intn = shape[ndim-1];
    CBLAS_INT incx = strides[ndim-1];
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        // The slices are vectors and `compute_slice_ptr` assumes matrices
        slice_ptr_a = compute_slice_ptr<T>(idx, A_data, ndim + 1, shape, strides);

        // XXX: introduce enum for fast paths & have general case as fallback?
        if (ord == 1) {
            res_data[idx] = _vector_norm_1<T, false>(slice_ptr_a, intn, incx);
        } else if (ord == -1) {
            res_data[idx] = _vector_norm_1<T, true>(slice_ptr_a, intn, incx);
        } else if (ord == 2) {
            res_data[idx] = _vector_norm_2<T, false>(slice_ptr_a, intn, incx);
        } else if (ord == -2) {
            res_data[idx] = _vector_norm_2<T, true>(slice_ptr_a, intn, incx);
        } else if (ord == std::numeric_limits<f64>::infinity()) {
            res_data[idx] = _matrix_norm_1<T, false>(slice_ptr_a, 1, intn, 0, incx);
        } else if (ord == -std::numeric_limits<f64>::infinity()) {
            res_data[idx] = _matrix_norm_1<T, true>(slice_ptr_a, 1, intn, 0, incx);
        } else {
            res_data[idx] = _vector_norm_general<T>(slice_ptr_a, ord, intn, incx);
        }
    }

    return {};
}


template <typename T> BatchStatus _matrix_norm(PyArrayObject *ap_Am, PyArrayObject *ap_res, Matrix_norm_ord ord) {
    using real_type = real_of_t<T>;

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

    bool needs_svd = (ord == Matrix_norm_ord::m_nuc || ord == Matrix_norm_ord::m_two || ord == Matrix_norm_ord::m_neg_two);

    // -------------------------------------------------------------------
    // Workspace allocation (when norms require SVD)
    // -------------------------------------------------------------------

    char jobu = 'N', jobvt = 'N';
    char range = (ord == Matrix_norm_ord::m_nuc) ? 'A' : 'I';
    CBLAS_INT m = shape[ndim - 2], n = shape[ndim - 1];
    CBLAS_INT min_mn = std::min(m, n);
    real_of_t<T> vl = 0, vu = 0;
    CBLAS_INT il = (ord == Matrix_norm_ord::m_two) ? 1 : min_mn;
    CBLAS_INT iu = il;
    CBLAS_INT ns = 0;
    T tmp = 0;
    CBLAS_INT info = 0;

    std::uint64_t rwork_size = 0;
    CBLAS_INT lwork = 0;
    if constexpr (is_complex_v<T>) { rwork_size = 17 * min_mn * min_mn; }
    Workspace<real_of_t<T>> rwork;
    Workspace<CBLAS_INT> iwork;

    if (needs_svd) {
        rwork = alloc_workspace<real_of_t<T>>(rwork_size);
        iwork = alloc_workspace<CBLAS_INT>(12 * min_mn);

        if constexpr (!is_complex_v<T>) {
            lapack::gesvdx(jobu, jobvt, range, m, n, NULL, m, vl, vu, il, iu, &ns, NULL, NULL, m, NULL, n, &tmp, -1, iwork.get(), &info);
        } else {
            lapack::gesvdx(jobu, jobvt, range, m, n, NULL, m, vl, vu, il, iu, &ns, NULL, NULL, m, NULL, n, &tmp, -1, rwork.get(), iwork.get(), &info);
        }

        if (info != 0) { return {BatchCode::lapack_error, "gesvdx", info}; }

        int64_t size = 0;
        lwork = _calc_lwork(tmp, 1.0, &size);
        if (lwork < 0) { return {BatchCode::too_large, "gesvdx", size}; }
    }

    /*
     * Allocate memory and chop buffer into parts
     *
     *   lwork   mn
     * |-------|----|
     * ^       ^
     * work    data_A
     */
    const std::uint64_t mn = (std::uint64_t)m * n;
    auto buf = alloc_workspace<T>(lwork + mn);
    auto S = alloc_workspace<real_of_t<T>>(min_mn + 1); // XXX: temporary shim, see gh-26824

    // Partition buffer
    T *cursor = buf.get();
    T *work = cursor; cursor += lwork;
    T *data_A = cursor;
    real_of_t<T> *data_S = S.get();


    // -------------------------------------------------------------------
    // Main batching loop
    // -------------------------------------------------------------------
    T *slice_ptr_A = NULL;
    real_of_t<T> res = 0;
    real_type *res_ptr = res_data;
    CBLAS_INT intm = shape[ndim - 2], intn = shape[ndim - 1];

    CBLAS_INT s0 = strides[ndim - 2];
    CBLAS_INT s1 = strides[ndim - 1];

    for (npy_intp idx = 0; idx < outer_size; idx++) {
        slice_ptr_A = compute_slice_ptr<T>(idx, A_data, ndim, shape, strides);

        switch (ord) {
            case Matrix_norm_ord::m_fro:
                res_ptr[idx] = _matrix_norm_F<T, false>(slice_ptr_A, intm, intn, s0, s1);
                break;

            case Matrix_norm_ord::m_nuc:
                copy_strided_to_f(data_A, slice_ptr_A, m, n, s0, s1, m);
                if constexpr (!is_complex_v<T>) {
                    lapack::gesvdx(jobu, jobvt, range, m, n, data_A, m, vl, vu, il, iu, &ns, data_S, NULL, m, NULL, n, work, lwork, iwork.get(), &info);
                } else {
                    lapack::gesvdx(jobu, jobvt, range, m, n, data_A, m, vl, vu, il, iu, &ns, data_S, NULL, m, NULL, n, work, lwork, rwork.get(), iwork.get(), &info);
                }
                res = 0;
                for (npy_intp i = 0; i < min_mn; i++) { res += data_S[i]; }
                res_ptr[idx] = res;
                break;

            case Matrix_norm_ord::m_one:
                res_ptr[idx] = _matrix_norm_1<T, false>(slice_ptr_A, intm, intn, s0, s1);
                break;

            case Matrix_norm_ord::m_neg_one:
                res_ptr[idx] = _matrix_norm_1<T, true>(slice_ptr_A, intm, intn, s0, s1);
                break;

            case Matrix_norm_ord::m_two:
                copy_strided_to_f(data_A, slice_ptr_A, m, n, s0, s1, m);
                if constexpr (!is_complex_v<T>) {
                    lapack::gesvdx(jobu, jobvt, range, m, n, data_A, m, vl, vu, il, iu, &ns, data_S, NULL, m, NULL, n, work, lwork, iwork.get(), &info);
                } else {
                    lapack::gesvdx(jobu, jobvt, range, m, n, data_A, m, vl, vu, il, iu, &ns, data_S, NULL, m, NULL, n, work, lwork, rwork.get(), iwork.get(), &info);
                }
                res_ptr[idx] = data_S[0];
                break;

            case Matrix_norm_ord::m_neg_two:
                copy_strided_to_f(data_A, slice_ptr_A, m, n, s0, s1, m);
                if constexpr (!is_complex_v<T>) {
                    lapack::gesvdx(jobu, jobvt, range, m, n, data_A, m, vl, vu, il, iu, &ns, data_S, NULL, m, NULL, n, work, lwork, iwork.get(), &info);
                } else {
                    lapack::gesvdx(jobu, jobvt, range, m, n, data_A, m, vl, vu, il, iu, &ns, data_S, NULL, m, NULL, n, work, lwork, rwork.get(), iwork.get(), &info);
                }
                res_ptr[idx] = data_S[0];
                break;

            case Matrix_norm_ord::m_inf:
                // Effectively transpose to iterate over rows instead of columns
                res_ptr[idx] = _matrix_norm_1<T, false>(slice_ptr_A, intn, intm, s1, s0);
                break;

            case Matrix_norm_ord::m_neg_inf:
                res_ptr[idx] = _matrix_norm_1<T, true>(slice_ptr_A, intn, intm, s1, s0);
                break;

            default: // Unknown enum
                return {BatchCode::lapack_error, "matrix_norm", -102};
        }
    }

    return {};
}



} // namespace sp_linalg
