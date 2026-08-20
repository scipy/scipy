#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False
"""
Part 2: directly uses scipy's blas_int without an intermediate .pxd layer.

This simulates a downstream package that directly adapts to SciPy's
blas_int type, using it throughout without int->blas_int conversion.
"""

from scipy.linalg.cython_blas cimport blas_int, ddot, daxpy, dgemm, dnrm2
from scipy.linalg.cython_lapack cimport dgetrf

from libc.stdlib cimport malloc, free


cpdef double direct_ddot(double[:] x, double[:] y):
    cdef:
        blas_int n = x.shape[0]
        blas_int incx = 1, incy = 1
    return ddot(&n, &x[0], &incx, &y[0], &incy)


cpdef direct_daxpy(double alpha, double[:] x, double[:] y):
    cdef:
        blas_int n = x.shape[0]
        blas_int incx = 1, incy = 1
    daxpy(&n, &alpha, &x[0], &incx, &y[0], &incy)


cpdef double direct_dnrm2(double[:] x):
    cdef:
        blas_int n = x.shape[0]
        blas_int incx = x.strides[0] // sizeof(double)
        double result = dnrm2(&n, &x[0], &incx)
    return result


cpdef direct_dgemm(double alpha, double[::1,:] a, double[::1,:] b,
                   double beta, double[::1,:] c):
    """Matrix multiply assuming Fortran-contiguous arrays."""
    cdef:
        char transa = b'N'
        char transb = b'N'
        blas_int m = a.shape[0]
        blas_int k = a.shape[1]
        blas_int n = b.shape[1]
        blas_int lda = m, ldb = k, ldc = m
    dgemm(&transa, &transb, &m, &n, &k,
          &alpha, &a[0, 0], &lda, &b[0, 0], &ldb, &beta, &c[0, 0], &ldc)


cpdef int direct_dgetrf(double[::1,:] a, int[:] ipiv_out):
    """LU factorization using blas_int directly for ipiv."""
    cdef:
        blas_int m = a.shape[0]
        blas_int n = a.shape[1]
        blas_int lda = m
        blas_int info
        blas_int *ipiv
        int min_mn = min(a.shape[0], a.shape[1])
        int i

    ipiv = <blas_int *>malloc(min_mn * sizeof(blas_int))
    if ipiv == NULL:
        return -1000

    dgetrf(&m, &n, &a[0, 0], &lda, ipiv, &info)

    for i in range(min_mn):
        ipiv_out[i] = <int>ipiv[i]

    free(ipiv)
    return <int>info


cpdef int direct_dgetrf_2(double[::1,:] a, blas_int[:] ipiv):
    """Similar to `direct_dgetrf` but accepts `ipiv` as a `blas_int` array."""
    cdef:
        blas_int m = a.shape[0]
        blas_int n = a.shape[1]
        blas_int lda = m
        blas_int info

    dgetrf(&m, &n, &a[0, 0], &lda, &ipiv[0], &info)
    return <int>info


cpdef int get_blas_int_size():
    """Return sizeof(blas_int) to verify correct type at compile+runtime."""
    return sizeof(blas_int)
