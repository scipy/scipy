"""
Routines for updating a Cholesky factorization

.. versionadded:: X.Y.Z

"""
#
# Copyright (C) 2022 Giacomo Meanti
# some functions copied from _decomp_update.pyx (Copyright 2014 Eric Moore)
# 
# Two algorithms are implemented in this module: one for rank-1 
# updates and downdates, which follows ref 1. and another for
# rank-p updates and downdates which follows ref 2.
#
# Useful references for Cholesky factorization up and downdating:
#
# 1. Pan, C. T. A modification to the linpack downdating algorithm. 
#    BIT Numerical Mathematics 30, 707-722 (1990).
# 2. Van de Geijn, R. A., Van Zee, F. G. High-Performance Up-and-Downdating 
#    via Householder-like Transformations. 
#    ACM Trans. Math. Softw. 38, 1-17 (2011).
#

__all__ = ["cholesky_update"]

cimport cython
from libc.math cimport sqrt, copysign, hypot
cimport libc.limits
cimport libc.float
cimport numpy as cnp
from . cimport cython_blas as blas_pointers
from . cimport cython_lapack as lapack_pointers

import numpy as np


# These are commented out in the numpy support we cimported above.
# Here I have declared them as taking void* instead of PyArrayDescr
# and object. In this file, only NULL is passed to these parameters.
cdef extern from *:
    cnp.ndarray PyArray_CheckFromAny(object, void*, int, int, int, void*)
    cnp.ndarray PyArray_FromArray(cnp.ndarray, void*, int)

# This is used in place of, e.g., cnp.NPY_C_CONTIGUOUS, to indicate that a C
# F or non contiguous array is acceptable.
DEF ARRAY_ANYORDER = 0


ctypedef float complex float_complex
ctypedef double complex double_complex
ctypedef fused float_t:
    float
    double
    float_complex
    double_complex

ctypedef fused real_t:
    float
    double


#------------------------------------------------------------------------------
# Helper functions for indexing and complex support
#------------------------------------------------------------------------------

cdef inline float_t* index2(float_t* a, int* astrides, int i, int j) nogil:
    return a + i*astrides[0] + j*astrides[1]

cdef inline float_t* index1(float_t* a, int* astrides, int i) nogil:
    return a + i*astrides[0]

cdef void float_t_conj(int n, float_t* x, int* xs) nogil:
    cdef int j
    if float_t is float_complex or float_t is double_complex:
        for j in range(n):
            index1(x, xs, j)[0] = index1(x, xs, j)[0].conjugate()

cdef void float_t_2d_conj(int m, int n, float_t* x, int* xs) nogil:
    cdef int i, j
    if float_t is float_complex or float_t is double_complex:
        for i in range(m):
            for j in range(n):
                index2(x, xs, i, j)[0] = index2(x, xs, i, j)[0].conjugate()

cdef float_t float_t_sqrt(float_t x) nogil:
    if float_t is float:
        return sqrt(x)
    elif float_t is double:
        return sqrt(x)
    elif float_t is float_complex:
        return <float_complex>sqrt(<double>((<float*>&x)[0]))
    else:
        return sqrt((<double*>&x)[0])

cdef inline bint float_t_less_than(float_t x, float_t y) nogil:
    if float_t is float or float_t is double:
        return x < y
    else:
        return x.real < y.real

#------------------------------------------------------------------------------
# BLAS and LAPACK wrappers
#------------------------------------------------------------------------------

cdef inline void scal(int n, float_t a, float_t* x, int incx) nogil:
    # Some? BLAS distributions don't properly handle negative stride.
    # https://stackoverflow.com/questions/44875236/why-are-non-positive-strides-disallowed-in-the-blas-gemm-family-of-functions
    if incx < 0:
        incx = -incx
        if float_t is float:
            blas_pointers.sscal(&n, &a, x - (n - 1) * incx, &incx)
        elif float_t is double:
            blas_pointers.dscal(&n, &a, x - (n - 1) * incx, &incx)
        elif float_t is float_complex:
            blas_pointers.cscal(&n, &a, x - (n - 1) * incx, &incx)
        else:  # float_t is double_complex:
            blas_pointers.zscal(&n, &a, x - (n - 1) * incx, &incx)
    else:
        if float_t is float:
            blas_pointers.sscal(&n, &a, x, &incx)
        elif float_t is double:
            blas_pointers.dscal(&n, &a, x, &incx)
        elif float_t is float_complex:
            blas_pointers.cscal(&n, &a, x, &incx)
        else:  # float_t is double_complex:
            blas_pointers.zscal(&n, &a, x, &incx)

cdef inline void axpy(int n, float_t a, float_t* x, int incx, float_t* y,
                      int incy) nogil:
    if float_t is float:
        blas_pointers.saxpy(&n, &a, x, &incx, y, &incy)
    elif float_t is double:
        blas_pointers.daxpy(&n, &a, x, &incx, y, &incy)
    elif float_t is float_complex:
        blas_pointers.caxpy(&n, &a, x, &incx, y, &incy)
    else:  # float_t is double_complex:
        blas_pointers.zaxpy(&n, &a, x, &incx, y, &incy)

cdef inline void larfg(int n, float_t* alpha, float_t* x, int incx,
                       float_t* tau) nogil:
    if float_t is float:
        lapack_pointers.slarfgp(&n, alpha, x, &incx, tau)
    elif float_t is double:
        lapack_pointers.dlarfgp(&n, alpha, x, &incx, tau)
    elif float_t is float_complex:
        lapack_pointers.clarfgp(&n, alpha, x, &incx, tau)
    else:  # float_t is double_complex:
        lapack_pointers.zlarfgp(&n, alpha, x, &incx, tau)

cdef inline void gemv(char* trans, int m, int n, float_t alpha, float_t* a,
                      int lda, float_t* x, int incx, float_t beta, float_t* y,
                      int incy) nogil:
    if float_t is float:
        blas_pointers.sgemv(trans, &m, &n, &alpha, a, &lda, x, &incx,
                &beta, y, &incy)
    elif float_t is double:
        blas_pointers.dgemv(trans, &m, &n, &alpha, a, &lda, x, &incx,
                &beta, y, &incy)
    elif float_t is float_complex:
        blas_pointers.cgemv(trans, &m, &n, &alpha, a, &lda, x, &incx,
                &beta, y, &incy)
    else:  # float_t is double_complex:
        blas_pointers.zgemv(trans, &m, &n, &alpha, a, &lda, x, &incx,
                &beta, y, &incy)

cdef inline void copy(int n, float_t* x, int incx, float_t* y, int incy) nogil:
    if float_t is float:
        blas_pointers.scopy(&n, x, &incx, y, &incy)
    elif float_t is double:
        blas_pointers.dcopy(&n, x, &incx, y, &incy)
    elif float_t is float_complex:
        blas_pointers.ccopy(&n, x, &incx, y, &incy)
    else:  # float_t is double_complex:
        blas_pointers.zcopy(&n, x, &incx, y, &incy)

cdef inline void ger(int m, int n, float_t alpha, float_t* x, int incx,
                     float_t* y, int incy, float_t* a, int lda) nogil:
    if float_t is float:
        blas_pointers.sger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda)
    elif float_t is double:
        blas_pointers.dger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda)
    elif float_t is float_complex:
        blas_pointers.cgerc(&m, &n, &alpha, x, &incx, y, &incy, a, &lda)
    else:
        blas_pointers.zgerc(&m, &n, &alpha, x, &incx, y, &incy, a, &lda)

cdef inline real_t nrm2(int n, real_t *x, int incx) nogil:
    if real_t is float:
        return blas_pointers.snrm2(&n, x, &incx)
    else:
        return blas_pointers.dnrm2(&n, x, &incx)

#------------------------------------------------------------------------------
# Negative larfg: generate a hausholder rotation with negative sign
#------------------------------------------------------------------------------

cdef inline real_t cathetus(real_t a, real_t b) nogil:
    return sqrt((a - b) * (a + b));

cdef void neg_larfg_real(int n, real_t *alpha, real_t *x, int incx,
                         real_t *tau) nogil:
    cdef real_t xnorm, beta

    xnorm = nrm2(n - 1, x, incx)
    beta = -copysign(cathetus(alpha[0], xnorm), alpha[0])
    tau[0] = (beta - alpha[0]) / beta
    scal(n - 1, 1 / (alpha[0] - beta), x, incx)
    alpha[0] = beta

cdef void neg_larfg_fc(int n, float_complex *alpha, float_complex *x, int incx,
                       float_complex *tau) nogil:
    cdef float xnorm, beta
    cdef int nm1 = n - 1
    xnorm = blas_pointers.scnrm2(&nm1, x, &incx)
    beta = -copysign(cathetus(alpha[0].real, xnorm), alpha[0].real)
    tau[0].real = (beta - alpha[0].real) / beta
    tau[0].imag = -alpha[0].imag / beta
    scal(n - 1, 1 / (alpha[0] - beta), x, incx)
    alpha[0] = beta

cdef void neg_larfg_dc(int n, double_complex *alpha, double_complex *x,
                       int incx, double_complex *tau) nogil:
    cdef double xnorm, beta
    cdef int nm1 = n - 1
    xnorm = blas_pointers.dznrm2(&nm1, x, &incx)
    beta = -copysign(cathetus(alpha[0].real, xnorm), alpha[0].real)
    tau[0].real = (beta - alpha[0].real) / beta
    tau[0].imag = -alpha[0].imag / beta
    scal(n - 1, 1 / (alpha[0] - beta), x, incx)
    alpha[0] = beta

cdef void neg_larfg(int n, float_t *alpha, float_t *x, int incx,
                    float_t *tau) nogil:
    if float_t is float or float_t is double:
        neg_larfg_real(n, alpha, x, incx, tau)
    elif float_t is float_complex:
        neg_larfg_fc(n, alpha, x, incx, tau)
    else:
        neg_larfg_dc(n, alpha, x, incx, tau)

#------------------------------------------------------------------------------
# Cholesky update (rank-1 and rank-p)
#------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef int upd1_u(float_t *r, int *rs, float_t *z, int *zs, int n,
                float_t eps, bint downdate) nogil:
    cdef float_t alpha = 1.0
    cdef float_t beta_prev = 1.0
    cdef float_t sign = -1.0 if downdate else 1.0
    cdef float_t a, beta
    cdef int i

    for i in range(n):
        # a = z[i] / r[i, i]
        a = index1(z, zs, i)[0] / index2(r, rs, i, i)[0]

        alpha = alpha + sign * (a.conjugate() * a)
        if float_t_less_than(alpha, eps):
            return -(i + 1)
        beta = float_t_sqrt(alpha)

        if n - i - 1 > 0:
            # z = -a * R + z
            float_t_conj(n - i - 1, index1(z, zs, i+1), zs)
            axpy(n=n - i - 1, a=-a.conjugate(), x=index2(r, rs, i, i + 1),
                 incx=rs[1], y=index1(z, zs, i + 1), incy=zs[0])

        scal(n - i, beta / beta_prev, index2(r, rs, i, i), rs[1])

        if n - i - 1 > 0:
            # R = -a/(beta*beta_m1) * z + R
            axpy(n=n - i - 1, a=sign * a / (beta * beta_prev),
                 x=index1(z, zs, i + 1), incx=zs[0],
                 y=index2(r, rs, i, i + 1), incy=rs[1])
            float_t_conj(n - i - 1, index1(z, zs, i+1), zs)
        beta_prev = beta

    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef int updp_u(float_t *r, int *rs, float_t *z, int *zs, float_t *w,
                int n, int p, bint zisF, float_t eps, bint downdate) nogil:
    """
    r is n x n square matrix, upper triangular
    z is n x p matrix, each column contains an update to r.
    w is n vector, contains working memory.

    https://www.cs.utexas.edu/users/flame/pubs/flawn41.pdf
    """
    cdef int i, j
    cdef float_t tau
    cdef char* T = 'T' if float_t is float or float_t is double else 'C'
    cdef char* N = 'N'
    cdef float_t sign = -1.0 if downdate else 1.0

    if zisF and zs[1] == 1 and p == 1:
        # doesn't matter as it's never actually used (since p == 1
        # goes to upd1_u), but otherwise lda parameter to ger is wrong
        zs[1] = n
    elif not zisF and zs[0] == 1 and p == 1:
        zs[0] = n

    for i in range(n):
        # Compute householder rotation. The rotation vector overwrites z[i, :],
        # and the scalar beta overwrites r[i, i].
        if downdate:
            # neg_larfg does not implement larfgp due to a loss of precision
            # in its naive implementation. Hence the signs are flipped and
            # we must flip all outputs at the end.
            neg_larfg(p + 1, index2(r, rs, i, i), x=index2(z, zs, i, 0),
                      incx=zs[1], tau=&tau)
            index2(r, rs, i, i)[0] = -index2(r, rs, i, i)[0]
            if tau != tau:
                return -(i + 1)
        else:
            larfg(p + 1, index2(r, rs, i, i), x=index2(z, zs, i, 0),
                  incx=zs[1], tau=&tau)
        if i >= n - 1:
            break

        # Working array names loosely following the paper:
        # r12: r[i, i+1:]   n - i - 1
        # c2: z[i+1:, :]    n - i - 1, p
        # u: z[i, :]        p
        # w:                n - i - 1

        # 1. copy r12 into w
        copy(n - i - 1, x=index2(r, rs, i, i+1), incx=rs[1], y=w, incy=1)
        # 2. w = tau * C2 @ u + tau * w
        if zisF:
            # w = tau * conj(C2) @ u + tau * w
            float_t_2d_conj(n - i - 1, p, index2(z, zs, i + 1, 0), zs)
            gemv(N, m=n-i-1, n=p, alpha=sign * tau, a=index2(z, zs, i+1, 0),
                 lda=zs[1], x=index2(z, zs, i, 0), incx=zs[1], beta=tau,
                 y=w, incy=1)
        else:
            # w = tau * C2^H @ u + tau * w
            gemv(T, m=p, n=n-i-1, alpha=sign * tau, a=index2(z, zs, i+1, 0),
                 lda=zs[0], x=index2(z, zs, i, 0), incx=zs[1], beta=tau,
                 y=w, incy=1)
        # 3. C2 = C2 - w @ u^H
        if zisF:
            # conj(C2) = conj(C2) - w @ u^H
            ger(m=n-i-1, n=p, alpha=-1.0, x=w, incx=1, y=index2(z, zs, i, 0),
                incy=zs[1], a=index2(z, zs, i+1, 0), lda=zs[1])
            float_t_2d_conj(n - i - 1, p, index2(z, zs, i + 1, 0), zs)
        else:
            # u @ w^H
            ger(m=p, n=n-i-1, alpha=-1.0, x=index2(z, zs, i, 0), incx=zs[1],
                y=w, incy=1, a=index2(z, zs, i+1, 0), lda=zs[0])
        # 4. r12 = r12 - w
        for j in range(n - i - 1):
            if downdate:  # refer to comment at start of loop.
                index2(r, rs, i, i + 1 + j)[0] = -index2(r, rs, i, i + 1 + j)[0] + w[j]
            else:
                index2(r, rs, i, i + 1 + j)[0] -= w[j]
    return 0


cdef validate_array(cnp.ndarray a, bint chkfinite):
    # here we check that a has positive strides and that its size is small
    # enough to fit in into an int, as BLAS/LAPACK require
    cdef bint copy = False
    cdef int j

    for j in range(a.ndim):
        if a.strides[j] <= 0:
            copy = True
        if (a.strides[j] / a.descr.itemsize) >= libc.limits.INT_MAX:
            copy = True
        if a.shape[j] >= libc.limits.INT_MAX:
            raise ValueError('Input array too large for use with BLAS')

    if chkfinite:
        if not np.isfinite(a).all():
            raise ValueError('array must not contain infs or NaNs')

    if copy:
        return PyArray_FromArray(a, NULL, cnp.NPY_F_CONTIGUOUS)
    return a


cdef tuple validate_inputs(object r0, object z0, bint overwrite_r,
                           bint overwrite_z,  int order_r, int order_z,
                           bint check_finite):
    cdef cnp.ndarray r, z
    cdef int typecode
    cdef int p

    order_r |= cnp.NPY_BEHAVED_NS | cnp.NPY_ELEMENTSTRIDES
    order_z |= cnp.NPY_BEHAVED_NS | cnp.NPY_ELEMENTSTRIDES
    if not overwrite_r:
        order_r |= cnp.NPY_ENSURECOPY
    if not overwrite_z:
        order_z |= cnp.NPY_ENSURECOPY

    # in the interests of giving better error messages take any number of
    # dimensions here.
    r = PyArray_CheckFromAny(r0, NULL, 0, 0, order_r, NULL)
    z = PyArray_CheckFromAny(z0, NULL, 0, 0, order_z, NULL)

    # Check shapes
    if z.ndim == 0:
        raise ValueError()

    if r.ndim != 2:
        raise ValueError()
    if r.shape[0] != r.shape[1]:
        raise ValueError()
    if r.shape[0] != z.shape[0]:
        raise ValueError()

    if z.ndim == 1:
        p = 1
    elif z.ndim == 2:
        p = z.shape[1]
    else:
        raise ValueError()

    # Check dtypes
    typecode = cnp.PyArray_TYPE(r)
    if typecode != cnp.PyArray_TYPE(z):
        raise ValueError('R and z must have the same dtype')

    if not (typecode == cnp.NPY_FLOAT or typecode == cnp.NPY_DOUBLE
            or typecode == cnp.NPY_CFLOAT or typecode == cnp.NPY_CDOUBLE):
        raise ValueError('Only arrays with dtypes float32, float64, '
                         'complex64, and complex128 are supported.')

    r = validate_array(r, check_finite)
    z = validate_array(z, check_finite)

    return r, z, typecode, r.shape[0], p


cdef void* extract(cnp.ndarray arr, int* arrs):
    with cython.cdivision(True):    # Assumes itemsize != 0.
        if arr.ndim == 2:
            arrs[0] = arr.strides[0] / cnp.PyArray_ITEMSIZE(arr)
            arrs[1] = arr.strides[1] / cnp.PyArray_ITEMSIZE(arr)
        elif arr.ndim == 1:
            arrs[0] = arr.strides[0] / cnp.PyArray_ITEMSIZE(arr)
            arrs[1] = 0
    return cnp.PyArray_DATA(arr)


@cython.embedsignature(True)
def cholesky_update(R, z, downdate=False, lower=False, overwrite_rz=False,
                    check_finite=True):
    """
    Rank-P Cholesky decomposition update.

    Given the Cholesky factorization :math:`A=LL^*` or :math:`A=U^*U`,
    returns the updated factorization of :math:`A + zz^*` or :math:`A - zz^*`.

    Parameters
    ----------
    R : (N, N) array_like
        Triangular matrix from the Cholesky decomposition of `A`.
    z : (N,) or (N, P) array_like
        Update vector(s)
    downdate : bool, optional
        Whether the rank-p update is to be added or subtracted from `A`.
    lower : bool, optional
        Whether `R` is a lower- or upper-triangular matrix.
    overwrite_rz : bool, optional
        Whether to overwrite data in `R` and `z`. If True, `R` will contain the
        updated factorization.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    R1 : (N, N) ndarray
        Updated Cholesky decomposition. Will have the same triangular pattern
        as `R`.

    See Also
    --------
    cholesky : Cholesky factorization

    Notes
    -----

    .. versionadded:: X.Y.Z

    References
    ----------
    .. [1] Van de Geijn, R. A., Van Zee, F. G. High-Performance
           Up-and-Downdating via Householder-like Transformations.
           ACM Trans. Math. Softw. 38, 1-17 (2011).
    .. [2] Pan, C. T. A modification to the linpack downdating algorithm. BIT
           Numerical Mathematics 30, 707-722 (1990).

    Raises
    ------
    LinAlgError : 
        If a downdate results in a non-positive definite matrix.

    Examples
    --------
    >>> from scipy import linalg
    >>> a = np.array([[  3., -1.,  0.],
    ...               [ -2.,  4.,  1.],
    ...               [ -3., -8.,  7.]])
    >>> r = linalg.cholesky(a, lower=False)

    Given this Cholesky decomposition, perform a rank-1 update

    >>> u = np.array([  1.,  4., -2.])
    >>> r_upd = linalg.cholesky_update(r, u, downdate=False)
    >>> r_upd
    array([[ 2.        ,  1.5       , -1.        ],
           [ 0.        ,  4.21307489, -1.30545982],
           [ 0.        ,  0.        ,  2.88023864]])

    The equivalent update can be obtained more slowly with the following.

    >>> a1 = a + np.outer(u, u)
    >>> r_upd_direct = linalg.cholesky(a1, lower=False)

    Check that the results are equivalent:

    >>> np.allclose(r_upd_direct, r_upd)
    True

    """
    cdef cnp.ndarray r1, z1
    cdef int typecode, n, p
    cdef bint chkfinite = check_finite, overwrite = overwrite_rz
    cdef bint is_upper = not lower, do_downdate = downdate

    cdef void* rptr
    cdef void* zptr
    cdef int rstrides[2]
    cdef int zstrides[2]

    r1, z1, typecode, n, p = validate_inputs(R, z, overwrite, overwrite,
        ARRAY_ANYORDER, ARRAY_ANYORDER, chkfinite)

    cdef float eps_f = n * libc.float.FLT_EPSILON
    cdef double eps_d = n * libc.float.DBL_EPSILON
    cdef float_complex eps_cf = eps_f
    cdef float_complex eps_cd = eps_d

    rptr = extract(r1, rstrides)
    cdef int ws
    if not is_upper:
        # Only upper is implemented. So for lower triangular inputs
        # we must transpose (just switch the memory layout since the
        # matrix is square).
        ws = rstrides[0]
        rstrides[0] = rstrides[1]
        rstrides[1] = ws

    cdef int info = 0
    cdef bint z_is_f
    cdef cnp.ndarray w0
    cdef void* wptr
    cdef cnp.npy_intp wlength = n
    if p == 1:  # z1 is 1D or column-vector
        zptr = extract(z1, zstrides)
        with nogil:
            if typecode == cnp.NPY_FLOAT:
                info = upd1_u(<float*>rptr, rstrides, <float*>zptr, zstrides,
                              n, eps_f, do_downdate)
            elif typecode == cnp.NPY_DOUBLE:
                info = upd1_u(<double*>rptr, rstrides, <double*>zptr, zstrides,
                              n, eps_d, do_downdate)
            elif typecode == cnp.NPY_CFLOAT:
                # TODO: Maybe only conj the triangle we're interested in?
                if not is_upper:
                    float_t_2d_conj(n, n, <float_complex*>rptr, rstrides)
                info = upd1_u(<float_complex*>rptr, rstrides,
                              <float_complex*>zptr, zstrides, n, eps_cf,
                              do_downdate)
                if not is_upper:
                    float_t_2d_conj(n, n, <float_complex*>rptr, rstrides)
            else:
                if not is_upper:
                    float_t_2d_conj(n, n, <double_complex*>rptr, rstrides)
                info = upd1_u(<double_complex*>rptr, rstrides,
                              <double_complex*>zptr, zstrides, n, eps_cd,
                              do_downdate)
                if not is_upper:
                    float_t_2d_conj(n, n, <double_complex*>rptr, rstrides)
    else:
        # Copy z if it's not contiguous
        if not cnp.PyArray_ISONESEGMENT(z1):
            z1 = PyArray_FromArray(z1, NULL, cnp.NPY_F_CONTIGUOUS)
        z_is_f = cnp.PyArray_IS_F_CONTIGUOUS(z1)
        zptr = extract(z1, zstrides)
        # Create work array (length n)
        w0 = cnp.PyArray_EMPTY(1, &wlength, typecode, 1)
        wptr = cnp.PyArray_DATA(w0)
        with nogil:
            if typecode == cnp.NPY_FLOAT:
                info = updp_u(<float*>rptr, rstrides, <float*>zptr, zstrides,
                              <float*>wptr, n, p, z_is_f, eps_f, do_downdate)
            elif typecode == cnp.NPY_DOUBLE:
                info = updp_u(<double*>rptr, rstrides, <double*>zptr, zstrides,
                              <double*>wptr, n, p, z_is_f, eps_f, do_downdate)
            elif typecode == cnp.NPY_CFLOAT:
                if not is_upper:
                    float_t_2d_conj(n, n, <float_complex*>rptr, rstrides)
                info = updp_u(<float_complex*>rptr, rstrides,
                              <float_complex*>zptr, zstrides,
                              <float_complex*>wptr, n, p, z_is_f, eps_f,
                              do_downdate)
                if not is_upper:
                    float_t_2d_conj(n, n, <float_complex*>rptr, rstrides)
            else:
                if not is_upper:
                    float_t_2d_conj(n, n, <double_complex*>rptr, rstrides)
                info = updp_u(<double_complex*>rptr, rstrides,
                              <double_complex*>zptr, zstrides,
                              <double_complex*>wptr, n, p, z_is_f, eps_f,
                              do_downdate)
                if not is_upper:
                    float_t_2d_conj(n, n, <double_complex*>rptr, rstrides)
    if info != 0:
        raise np.linalg.LinAlgError(
            "Update causes %d-th leading minor of the "
            "array to not be positive definite." % (-info))
    return r1

cnp.import_array()