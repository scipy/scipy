#pragma once

#include "Python.h"
#include <limits>
#include <climits>
#include <cassert>
#include "numpy/arrayobject.h"
#include "npy_cblas.h"
#include "_lapack_trampolines.h"
#include "_npymath.h"


using namespace _numpymath;


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
    char *data_ptr;

    npy_intp outer_size;  // math.prod(a.shape[:-2])
    npy_intp m, n;        // core dimensions

    iter_data_t(PyArrayObject *arr) {
        ndim = PyArray_NDIM(arr);
        shape = PyArray_SHAPE(arr);
        strides = PyArray_STRIDES(arr);
        data_ptr = (char *)PyArray_DATA(arr);

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
};


/*
 * Copy the slice number `idx` into the `buffer`.
 */
template<typename T>
void copy_slice(const iter_data_t *iter_data, const npy_intp idx, T *buffer) {

    assert((0 <= idx) && (idx < iter_data->outer_size));

    // compute the offset for the current slice
    npy_intp offset = 0;
    npy_intp temp_idx = idx;
    for (npy_intp i = iter_data->ndim - 3; i >= 0; i--) {
        offset += (temp_idx % iter_data->shape[i]) * iter_data->strides[i];
        temp_idx /= iter_data->shape[i];
    }

    char *slice_ptr = iter_data->data_ptr + offset;

    // core shape is (m, n)
    npy_intp m = iter_data->m;
    npy_intp n = iter_data->n;

    // core strides are (s2, s1)
    npy_intp s2 = iter_data->strides[iter_data->ndim - 2];
    npy_intp s1 = iter_data->strides[iter_data->ndim - 1];

    // copy the current m-by-n slice into the temp storage *in the Fortran order*
    for(npy_intp i=0; i < m; i++) {
        for(npy_intp j=0; j < n; j++) {
            buffer[i + j*iter_data->m] = *(T *)(slice_ptr + i*s2 + j*s1);
        }
    }
}



/***********************************
 *   Inversion via GETRF / GETRI
 ***********************************/

/*
 * Invert a 2D slice.
 *
 * Return is 0 on success; on failure, return the failing operation's
 * `info` variable.
 *
 */
template<typename T>
int
invert_slice(CBLAS_INT N, T *A, CBLAS_INT lda, CBLAS_INT *ipiv, T *work, CBLAS_INT lwork) {
    CBLAS_INT info;

    // factorize
    GETRF(&N, &N, A, &lda, ipiv, &info);
    if (info != 0) {
        return info;
    }

    // compute the inverse
    GETRI(&N, A, &lda, ipiv, work, &lwork, &info);
    return info;
}


/*
 * The main batched-inv routine, to be called from the python API layer.
 *
 * This is a computational routine, and it assumes that the caller has done all input
 * validation before calling it.
 * In particular:
 *
 * For the input array `arr`:
 *
 * - The caller must check that the dtype is LAPACK-compatible before calling this routine.
 * - There are no restriction on the input array strides, this routine will handle
 *   them correctly.
 *
 * For the output array, `arr_inv`, there are two options:
 * 
 *   1. It is a C-ordered array of the same shape and dtype as `arr`. Then it is filled
 *      with the inverse of `arr`.
 *   2. It coincides with `arr`. Then the inversion will be made in-place. 
 *
 * Currently, the latter behavior (a.k.a. `overwrite_a=True`) is only available for
 * 2D arrays. Again, the caller is responsible for checking that it is F-ordered.
 *
 * Return value is:
 *  - 0 on success;
 *  - LLONG_MIN on a memory error;
 *  - LAPACK info variable on a LAPACK error
 */
template<typename T>
inline int inv_loop(PyArrayObject *arr, PyArrayObject *arr_inv)
{
    /*
     * Check if we're reusing the memory. The caller has already checked that
     * the `arr` array has compatible strides, so we are not checking it here again.
     * In particular, the core dimensions *must* be F-ordered.
     *
     */
    bool overwrite_a = (arr == arr_inv);

    // looping/output variables
    bool all_failed=true;
    long long status=0, is_ok;
    T *ret_data = (T *)PyArray_DATA(arr_inv);

    /*
     * Prepare the data for looping over the batch dimensions.
     */
    iter_data_t iter_data(arr);
    npy_intp n = iter_data.n;    // core dimensions are (n, n)

    /* 
     * Workspace query. GETRI calculates this by passing LWORK = -1.
     * As query is meant to get optimal work size and not for actual decomposition,
     * no need to pass matrices.
     */
    CBLAS_INT N=(CBLAS_INT)n, lda=N, lwork=-1;
    CBLAS_INT info;
    T wrk;

    GETRI(&N, NULL, &lda, NULL, &wrk, &lwork, &info);

    /*
     * The factor of 1.01 here mirrors 
     * https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/_basic.py#L1154
     *
     * It was added in commit
     * https://github.com/scipy/scipy/commit/dfb543c147c
     * to avoid a "curious segfault with 500x500 matrices and OpenBLAS".
     */
    lwork = (CBLAS_INT)(1.01 * real_part(wrk));

    // Allocate work and other arrays.
    T *A;
    if (overwrite_a) {
        A = (T *)PyArray_DATA(arr_inv);
    } else {
        A = (T *)malloc(n*n*sizeof(T));
    }
    CBLAS_INT *ipiv = (CBLAS_INT *)malloc(n*sizeof(CBLAS_INT));  
    T *work = (T *)malloc(lwork*sizeof(T));

    if ((A == NULL) || (ipiv == NULL) || (work == NULL)) {
        PyErr_NoMemory();
        status = LLONG_MIN;
        goto done;
    }

    /* Finally, proceed to inverting the input matrix */

    /*
     * `overwrite_a=True` : 2D only, no looping, no copying.
     */
    if(overwrite_a){
        status = invert_slice(N, A, lda, ipiv, work, lwork);
        goto done;
    }

    /*
     * Main NxN slice loop
     */
    for(npy_intp idx=0; idx < iter_data.outer_size; idx++) {

        // fill the buffer with the current slice
        copy_slice(&iter_data, idx, A);

        // call LAPACK on the current slice
        is_ok = invert_slice(N, A, lda, ipiv, work, lwork);

        all_failed &= (is_ok != 0);
        if(is_ok == 0) {
            // copy to the output buffer, swap CF
            swap_cf(A, ret_data + idx*n*n, n, n, n);
        }
        else {
            // Either GETRF or GETRI failed.
            nan_matrix(ret_data + idx*n*n, n, n);
            status = is_ok;
        }
    }

    // XXX make it talk to np.errstate?

 done:
    if (!overwrite_a) {
        free(A);
    }
    free(ipiv);
    free(work);
    return all_failed ? status : 0;

};

