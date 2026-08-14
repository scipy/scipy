#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False
"""
Part 1: sklearn-like BLAS/LAPACK wrappers with int->blas_int conversion.

This module includes _blas_int.pxi (generated from the .pxi.in template)
which defines blas_int as either `int` (LP64) or `int64_t` (ILP64).
All wrapper functions accept `int` parameters and convert to `blas_int`
before calling the actual BLAS/LAPACK routines.
"""

include "_blas_int.pxi"

from scipy.linalg.cython_blas cimport ddot, daxpy, dgemm, dnrm2
from scipy.linalg.cython_lapack cimport dgetrf

from libc.stdlib cimport malloc, free

# Pointer casts from `const double *` to `double *` below are needed because
# SciPy's cython_blas.pxd lacks const qualifiers. This pattern is copied from
# scikit-learn's sklearn/utils/_cython_blas.pyx.
# See: https://github.com/scipy/scipy/issues/14262


cdef double _ddot(int n, const double *x, int incx,
                  const double *y, int incy) noexcept nogil:
    cdef blas_int bn = n, bincx = incx, bincy = incy
    return ddot(&bn, <double *>x, &bincx, <double *>y, &bincy)


cdef void _daxpy(int n, double alpha, const double *x, int incx,
                 double *y, int incy) noexcept nogil:
    cdef blas_int bn = n, bincx = incx, bincy = incy
    daxpy(&bn, &alpha, <double *>x, &bincx, y, &bincy)


cdef double _dnrm2(int n, const double *x, int incx) noexcept nogil:
    cdef:
        blas_int bn = n
        blas_int bincx = incx
        double nrm2 = dnrm2(&bn, <double *>x, &bincx)
    return nrm2


cdef void _dgemm(char *transa, char *transb, int m, int n, int k,
                 double alpha, const double *a, int lda,
                 const double *b, int ldb,
                 double beta, double *c, int ldc) noexcept nogil:
    cdef blas_int bm = m, bn = n, bk = k
    cdef blas_int blda = lda, bldb = ldb, bldc = ldc
    dgemm(transa, transb, &bm, &bn, &bk,
          &alpha, <double *>a, &blda, <double *>b, &bldb, &beta, c, &bldc)


cdef int _dgetrf(int m, int n, double *a, int lda,
                 int *ipiv) noexcept nogil:
    """LU factorization. The tricky part: ipiv is an output array of
    blas_int, but the downstream API uses int. We must allocate a
    temporary blas_int array, call LAPACK, then copy back."""
    cdef:
        blas_int bm = m, bn = n, blda = lda, info
        blas_int *bipiv
        int min_mn = min(m, n)
        int i

    bipiv = <blas_int *>malloc(min_mn * sizeof(blas_int))
    if bipiv == NULL:
        return -1000  # allocation failure

    dgetrf(&bm, &bn, a, &blda, bipiv, &info)

    for i in range(min_mn):
        ipiv[i] = <int>bipiv[i]

    free(bipiv)
    return <int>info


cpdef int get_blas_int_size():
    """Return sizeof(blas_int) to verify correct type at compile+runtime."""
    return sizeof(blas_int)
