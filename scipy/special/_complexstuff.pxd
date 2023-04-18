# -*-cython-*-
#
# Common functions required when doing complex arithmetic with Cython.
#

import cython
cimport numpy as np
from libc.math cimport (
    isnan, isinf, isfinite, log1p, fabs, exp, log, sin, cos, sqrt, pow
)

cdef extern from "_complexstuff.h":
    double npy_cabs(np.npy_cdouble z) nogil
    double npy_carg(np.npy_cdouble z) nogil
    np.npy_cdouble npy_clog(np.npy_cdouble z) nogil
    np.npy_cdouble npy_cexp(np.npy_cdouble z) nogil
    np.npy_cdouble npy_csin(np.npy_cdouble z) nogil
    np.npy_cdouble npy_ccos(np.npy_cdouble z) nogil
    np.npy_cdouble npy_csqrt(np.npy_cdouble z) nogil
    np.npy_cdouble npy_cpow(np.npy_cdouble x, np.npy_cdouble y) nogil

DEF tol = 2.220446092504131e-16

ctypedef double complex double_complex

ctypedef fused number_t:
    double
    double_complex

ctypedef union _complex_pun:
    np.npy_cdouble npy
    double_complex c99

cdef inline np.npy_cdouble npy_cdouble_from_double_complex(
        double_complex x) noexcept nogil:
    cdef _complex_pun z
    z.c99 = x
    return z.npy

cdef inline double_complex double_complex_from_npy_cdouble(
        np.npy_cdouble x) noexcept nogil:
    cdef _complex_pun z
    z.npy = x
    return z.c99

cdef inline bint zisnan(number_t x) noexcept nogil:
    if number_t is double_complex:
        return isnan(x.real) or isnan(x.imag)
    else:
        return isnan(x)

cdef inline bint zisfinite(number_t x) noexcept nogil:
    if number_t is double_complex:
        return isfinite(x.real) and isfinite(x.imag)
    else:
        return isfinite(x)

cdef inline bint zisinf(number_t x) noexcept nogil:
    if number_t is double_complex:
        return not zisnan(x) and not zisfinite(x)
    else:
        return isinf(x)

cdef inline double zreal(number_t x) noexcept nogil:
    if number_t is double_complex:
        return x.real
    else:
        return x

cdef inline double zabs(number_t x) noexcept nogil:
    if number_t is double_complex:
        return npy_cabs(npy_cdouble_from_double_complex(x))
    else:
        return fabs(x)

cdef inline double zarg(double complex x) noexcept nogil:
    return npy_carg(npy_cdouble_from_double_complex(x))

cdef inline number_t zlog(number_t x) noexcept nogil:
    cdef np.npy_cdouble r
    if number_t is double_complex:
        r = npy_clog(npy_cdouble_from_double_complex(x))
        return double_complex_from_npy_cdouble(r)
    else:
        return log(x)

cdef inline number_t zexp(number_t x) noexcept nogil:
    cdef np.npy_cdouble r
    if number_t is double_complex:
        r = npy_cexp(npy_cdouble_from_double_complex(x))
        return double_complex_from_npy_cdouble(r)
    else:
        return exp(x)

cdef inline number_t zsin(number_t x) noexcept nogil:
    cdef np.npy_cdouble r
    if number_t is double_complex:
        r = npy_csin(npy_cdouble_from_double_complex(x))
        return double_complex_from_npy_cdouble(r)
    else:
        return sin(x)

cdef inline number_t zcos(number_t x) noexcept nogil:
    cdef np.npy_cdouble r
    if number_t is double_complex:
        r = npy_ccos(npy_cdouble_from_double_complex(x))
        return double_complex_from_npy_cdouble(r)
    else:
        return cos(x)

cdef inline number_t zsqrt(number_t x) noexcept nogil:
    cdef np.npy_cdouble r
    if number_t is double_complex:
        r = npy_csqrt(npy_cdouble_from_double_complex(x))
        return double_complex_from_npy_cdouble(r)
    else:
        return sqrt(x)

cdef inline number_t zpow(number_t x, double y) noexcept nogil:
    cdef np.npy_cdouble r, z
    # FIXME
    if number_t is double_complex:
        z.real = y
        z.imag = 0.0
        r = npy_cpow(npy_cdouble_from_double_complex(x), z)
        return double_complex_from_npy_cdouble(r)
    else:
        return pow(x, y)

cdef inline double_complex zpack(double zr, double zi) noexcept nogil:
    cdef np.npy_cdouble z
    z.real = zr
    z.imag = zi
    return double_complex_from_npy_cdouble(z)

@cython.cdivision(True)
cdef inline double complex zlog1(double complex z) noexcept nogil:
    """
    Compute log, paying special attention to accuracy around 1. We
    implement this ourselves because some systems (most notably the
    Travis CI machines) are weak in this regime.

    """
    cdef:
        int n
        double complex coeff = -1
        double complex res = 0

    if zabs(z - 1) > 0.1:
        return zlog(z)
    z = z - 1
    if z == 0:
        return 0
    for n in range(1, 17):
        coeff *= -z
        res += coeff/n
        if zabs(res/coeff) < tol:
            break
    return res
