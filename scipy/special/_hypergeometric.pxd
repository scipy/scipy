from libc.math cimport fabs, exp, floor, M_PI

from . cimport sf_error
from ._cephes cimport poch

cdef extern from "numpy/npy_math.h":
    double npy_isnan(double x) nogil
    double NPY_NAN
    double NPY_INFINITY

cdef extern from 'specfun_wrappers.h':
    double hypU_wrap(double, double, double) nogil
    double hyp1f1_wrap(double, double, double) nogil

DEF EPS = 2.220446049250313e-16


cdef inline double hyperu(double a, double b, double x) nogil:
    if x < 0.0:
        sf_error.error("hyperu", sf_error.DOMAIN, NULL)
        return NPY_NAN

    if x == 0.0:
        if b > 1.0:
            # DMLF 13.2.16-18
            sf_error.error("hyperu", sf_error.SINGULAR, NULL)
            return NPY_INFINITY
        else:
            # DLMF 13.2.14-15 and 13.2.19-21
            return poch(1.0 - b + a, -a)

    return hypU_wrap(a, b, x)


cdef inline double hyp1f1(double a, double b, double x) nogil:
    if npy_isnan(a) or npy_isnan(b) or npy_isnan(x):
        return NPY_NAN
    if b <= 0 and b == floor(b):
        return NPY_INFINITY
    elif a == 0 or x == 0:
        return 1
    elif a == -1:
        return 1 - x / b
    elif a == b:
        return exp(x)
    elif a - b == 1:
        return (1 + x / b) * exp(x)
    elif a == 1 and b == 2:
        return (exp(x) - 1) / x
    elif a <= 0 and a == floor(a):
        # For `a` a negative integer the series is finite.
        return hyp1f1_series(a, b, x)

    if b > 0 and (fabs(a) + 1) * fabs(x) < 0.9 * b:
        # At for the kth term of the series we are multiplying by
        #
        # t_k = (a + b) * x / ((b + k) * (k + 1))
        #
        # We have that
        #
        # |t_k| < (|a| + 1) * |x| / |b|,
        #
        # which means that in this branch we get geometric
        # convergence.
        return hyp1f1_series(a, b, x)

    return hyp1f1_wrap(a, b, x)


cdef inline double hyp1f1_series(double a, double b, double x) nogil:
    cdef int k
    cdef double term = 1
    cdef double result = 1
    for k in range(500):
        term *= (a + k) * x / (b + k) / (k + 1)
        result += term
        if fabs(term) <= EPS * fabs(result):
            break
    else:
        sf_error.error("hyp1f1", sf_error.NO_RESULT, NULL)
        result = NPY_NAN
    return result
