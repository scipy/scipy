#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False
"""
Part 1 consumer: cimports from the internal _blas_lapack_wrappers module.

This validates that the internal .pxd works correctly for cross-module
cimport, which is the pattern used by scikit-learn (e.g., other Cython
modules cimporting from sklearn.utils._cython_blas).
"""

from libc.limits cimport INT_MAX

from ilp64_test_package._blas_lapack_wrappers cimport (
    _ddot, _daxpy, _dnrm2, _dgemm, _dgetrf, get_blas_int_size
)


cpdef double consumer_ddot(double[:] x, double[:] y):
    if x.shape[0] > INT_MAX:
        raise ValueError("Integer overflow in consumer_ddot.")
    cdef int n = x.shape[0]
    return _ddot(n, &x[0], 1, &y[0], 1)


# Strictly speaking, all functions should check for integer overflow of the input array
# size. We are omitting these checks below, for brevity.


cpdef consumer_daxpy(double alpha, double[:] x, double[:] y):
    cdef int n = x.shape[0]
    _daxpy(n, alpha, &x[0], 1, &y[0], 1)


cpdef consumer_dnrm2(double[:] x):
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(double)
        double nrm2 = _dnrm2(n, &x[0], incx)
    return nrm2


cpdef consumer_dgemm(double alpha, double[::1,:] a, double[::1,:] b,
                     double beta, double[::1,:] c):
    """Matrix multiply assuming Fortran-contiguous arrays."""
    cdef:
        char transa = b'N'
        char transb = b'N'
        int m = a.shape[0]
        int k = a.shape[1]
        int n = b.shape[1]
    _dgemm(&transa, &transb, m, n, k,
           alpha, &a[0, 0], m, &b[0, 0], k, beta, &c[0, 0], m)


cpdef int consumer_dgetrf(double[::1,:] a, int[:] ipiv):
    """LU factorization via the internal wrappers."""
    cdef int m = a.shape[0], n = a.shape[1]
    return _dgetrf(m, n, &a[0, 0], m, &ipiv[0])


cpdef int consumer_blas_int_size():
    return get_blas_int_size()
