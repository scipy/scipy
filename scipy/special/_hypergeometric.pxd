from libc.math cimport fabs, exp, floor, M_PI

import cython

from . cimport sf_error
from ._cephes cimport expm1, poch

cdef extern from "numpy/npy_math.h":
    double npy_isnan(double x) nogil
    double NPY_NAN
    double NPY_INFINITY

cdef extern from 'specfun_wrappers.h':
    double hypU_wrap(double, double, double) nogil
    double hyp1f1_wrap(double, double, double) nogil

DEF EPS = 2.220446049250313e-16
DEF ACCEPTABLE_RTOL = 1e-7


@cython.cdivision(True)
cdef inline double hyperu(double a, double b, double x) nogil:
    if npy_isnan(a) or npy_isnan(b) or npy_isnan(x):
        return NPY_NAN

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


@cython.cdivision(True)
cdef inline double hyp1f1(double a, double b, double x) nogil:
    if npy_isnan(a) or npy_isnan(b) or npy_isnan(x):
        return NPY_NAN
    if b <= 0 and b == floor(b):
        # There is potentially a pole.
        if b <= a < 0 and a == floor(a):
            # The Pochammer symbol (a)_n cancels the pole.
            return hyp1f1_series_track_convergence(a, b, x)
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
        return expm1(x) / x
    elif a <= 0 and a == floor(a):
        # The geometric series is finite in this case, but it could
        # still suffer from cancellation.
        return hyp1f1_series_track_convergence(a, b, x)

    if b > 0 and (fabs(a) + 1) * fabs(x) < 0.9 * b:
        # For the kth term of the series we are multiplying by
        #
        # t_k = (a + k) * x / ((b + k) * (k + 1))
        #
        # We have that
        #
        # |t_k| < (|a| + 1) * |x| / |b|,
        #
        # which means that in this branch we get geometric
        # convergence.
        return hyp1f1_series(a, b, x)

    return hyp1f1_wrap(a, b, x)


@cython.cdivision(True)
cdef inline double hyp1f1_series_track_convergence(
    double a,
    double b,
    double x
) nogil:
    # The hypergeometric series can suffer from cancelation or take a
    # prohibitive number of terms to converge. This function computes
    # the series while monitoring those conditions.
    cdef int k
    cdef double apk, bpk
    cdef double term = 1
    cdef double result = 1
    cdef double abssum = result
    for k in range(1000):
        apk = a + k
        bpk = b + k
        if bpk != 0:
            term *= apk * x / bpk / (k + 1)
        elif apk == 0:
            # The Pochammer symbol in the denominator has become zero,
            # but we still have the continuation formula DLMF 13.2.5.
            term = 0
        else:
            # We hit a pole
            return NPY_NAN
        abssum += fabs(term)
        result += term
        if fabs(term) <= EPS * fabs(result):
            break
    else:
        sf_error.error("hyp1f1", sf_error.NO_RESULT, NULL)
        return NPY_NAN

    if k * EPS * abssum <= ACCEPTABLE_RTOL * fabs(result):
        return result
    sf_error.error("hyp1f1", sf_error.NO_RESULT, NULL)
    return NPY_NAN


@cython.cdivision(True)
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
