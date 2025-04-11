#pragma once

#include "Python.h"
#include <limits>
#include <cassert>
#include "numpy/arrayobject.h"
#include "npy_cblas.h"
#include "_lapack_trampolines.h"


template<typename T> struct numeric_limits {};

template<>
struct numeric_limits<float>{
    static constexpr double one = 1.0f;
    static constexpr float nan = std::numeric_limits<float>::quiet_NaN();
};

template<>
struct numeric_limits<double>{
    static constexpr double one = 1.0;
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
};


template<>
struct numeric_limits<npy_cfloat>{
    static constexpr npy_cfloat one = {1.0f, 0.0f};
    static constexpr npy_cfloat nan = {std::numeric_limits<float>::quiet_NaN(),
                                       std::numeric_limits<float>::quiet_NaN()};
};

template<>
struct numeric_limits<npy_cdouble>{
    static constexpr npy_cdouble one = {1.0, 0.0};
    static constexpr npy_cdouble nan = {std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN()};
};


/* 
 * Helpers for filling/rearranging matrices
 */


/* identity square matrix generation */
template<typename typ>
static inline void
identity_matrix(typ *matrix, Py_ssize_t n)
{
    Py_ssize_t i;
    /* in IEEE floating point, zeroes are represented as bitwise 0 */
    memset((void *)matrix, 0, n*n*sizeof(typ));

    for (i = 0; i < n; ++i)
    {
        *matrix = numeric_limits<typ>::one;
        matrix += n+1;
    }
}


/* m-by-n matrix full of nans */
template<typename typ>
static inline void
nan_matrix(typ *matrix, Py_ssize_t m, Py_ssize_t n)
{
    for (Py_ssize_t i=0; i < m*n; i++) {
        *(matrix + i) = numeric_limits<typ>::nan;
    }
}


// parroted from sqrtm
template<typename T>
inline void
swap_cf(T* src, T* dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    T *bb = dst;
    T *aa = src;
    if (r < 16) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf(src, dst, r2, c, n);
            swap_cf(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf(src, dst, r, c2, n);
            swap_cf(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
}


/*
 * Visit (m, n)-shaped slices of a ndim>=3 array in C order.
 * Copy each slice into a buffer, make the buffer F-ordered
 * for LAPACK.
 */
struct iter_data_t
{
    // grab useful data from the array to iterate
    npy_intp ndim;
    npy_intp *shape;
    npy_intp *strides;
    void *data_ptr;

    npy_intp outer_size;  // math.prod(a.shape[:-2])
    npy_intp m, n;        // core dimensions

    iter_data_t(PyArrayObject *arr) {
        ndim = PyArray_NDIM(arr);
        shape = PyArray_SHAPE(arr);
        strides = PyArray_STRIDES(arr);
        data_ptr = PyArray_DATA(arr);

        m = shape[ndim - 2];
        n = shape[ndim - 1];

        /* outer_size: the number of slices in the batch, math.prod(a.shape[:-2] */
        outer_size = 1;
        if(ndim > 2) {
            for(int i=0; i < ndim - 2; i++) {
                outer_size *= shape[i];
            }
        }
    }

    template<typename T>
    void get_slice(const npy_intp idx, T *buffer) {

        assert((0 <= idx) && (idx < outer_size));

        /* parroted from sqrtm*/
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }

        T *slice_ptr = (T *)(data_ptr) + offset/sizeof(T);

        // copy the current n-x-n slice into the temp storage *in the Fortran order*
        for(npy_intp i=0; i < m; i++) {
            for(npy_intp j=0; j < n; j++) {
                buffer[i + j*m] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(T)) + (j*strides[ndim - 1]/sizeof(T)));
            }
        }
    }
};



template<typename T>
inline void inv_loop2(PyArrayObject *a, PyArrayObject *a_inv)
{
    /*
     * Prepate the data for looping over the batch dimensions.
     */
    iter_data_t iter_data(a);
    npy_intp n = iter_data.n;    // core dimensions are (n, n)

    gesv_data_t<T> gesv_data((fortran_int)n);
    T *ret_data = (T *)PyArray_DATA(a_inv);

    /*
     * Main NxN slice loop
     */
    for(npy_intp idx=0; idx < iter_data.outer_size; idx++) {
        T *buffer = gesv_data.a;
        iter_data.get_slice(idx, buffer);   // fill the buffer with the current slice

        /* Call LAPACK  */
        identity_matrix(gesv_data.b, n);
        call_gesv(gesv_data);

        if (gesv_data.info != 0) {
            // XXX raise or quietly fill with nans
            for(npy_intp i=0; i < n*n; i++) {
                gesv_data.b[i] = std::numeric_limits<T>::quiet_NaN();
            }
        }

        // copy to the output buffer, swap CF; XXX: can use dcopy from BLAS?
        swap_cf(gesv_data.b, ret_data + idx*n*n, n, n, n);
    }
};



template<typename T>
inline int inv_loop(PyArrayObject *a, PyArrayObject *a_inv)
{
    /*
     * Prepate the data for looping over the batch dimensions.
     */
    iter_data_t iter_data(a);
    npy_intp n = iter_data.n;    // core dimensions are (n, n)

    getrf_data_t<T> getrf_data((fortran_int)n, (fortran_int)n);
    if(getrf_data.info != 0) {
        return getrf_data.info;
    }

    getri_data_t<T> getri_data(getrf_data);
    if(getri_data.info != 0) {
        return getri_data.info;
    }

    T *ret_data = (T *)PyArray_DATA(a_inv);

    /*
     * Main NxN slice loop
     */
    for(npy_intp idx=0; idx < iter_data.outer_size; idx++) {
        T *buffer = getrf_data.a;
        iter_data.get_slice(idx, buffer);   // fill the buffer with the current slice

        // factorize the current slice 
        call_getrf(getrf_data);

        if (getrf_data.info != 0) {
            // LAPACK error.
            // XXX raise or quietly fill with nans
            std::cout << "problem at idx=" << idx << "\n";

            nan_matrix(ret_data + idx*n*n, getri_data.n, getri_data.n);
            continue; 
        }

        // prepare the data for the GETRI call
        // (other getri_data members do not change between iterations
        getri_data.a = getrf_data.a;
        getri_data.ipiv = getrf_data.ipiv;
        call_getri(getri_data);

        if (getri_data.info != 0) {
            // LAPACK error.
            // XXX raise or quietly fill with nans
            std::cout << "GETRI info = " << getri_data.info << "\n";

            nan_matrix(ret_data + idx*n*n, getri_data.n, getri_data.n);
            continue;
        }

        // copy to the output buffer, swap CF; XXX: can use dcopy from BLAS?
        swap_cf(getri_data.a, ret_data + idx*n*n, n, n, n);
    }
    return 0;
};

