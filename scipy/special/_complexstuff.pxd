# -*-cython-*-
#
# Common functions required when doing complex arithmetic with Cython.
#

cimport numpy as np
cimport libc.math

cdef extern from "_complexstuff.h":
    double npy_cabs(np.npy_cdouble z) nogil
    np.npy_cdouble npy_clog(np.npy_cdouble z) nogil
    np.npy_cdouble npy_cexp(np.npy_cdouble z) nogil
    double npy_log1p(double x) nogil
    int npy_isnan(double x) nogil
    int npy_isinf(double x) nogil
    int npy_isfinite(double x) nogil
    double inf "NPY_INFINITY"
    double pi "NPY_PI"
    double nan "NPY_NAN"

ctypedef double complex double_complex

ctypedef fused number_t:
    double
    double_complex

ctypedef union _complex_pun:
    np.npy_cdouble npy
    double_complex c99

cdef inline np.npy_cdouble npy_cdouble_from_double_complex(
        double_complex x) nogil:
    cdef _complex_pun z
    z.c99 = x
    return z.npy

cdef inline double_complex double_complex_from_npy_cdouble(
        np.npy_cdouble x) nogil:
    cdef _complex_pun z
    z.npy = x
    return z.c99

cdef inline bint zisnan(number_t x) nogil:
    if number_t is double_complex:
        return npy_isnan(x.real) or npy_isnan(x.imag)
    else:
        return npy_isnan(x)

cdef inline bint zisfinite(number_t x) nogil:
    if number_t is double_complex:
        return npy_isfinite(x.real) and npy_isfinite(x.imag)
    else:
        return npy_isfinite(x)

cdef inline bint zisinf(number_t x) nogil:
    if number_t is double_complex:
        return not zisnan(x) and not zisfinite(x)
    else:
        return npy_isinf(x)

cdef inline double zabs(number_t x) nogil:
    if number_t is double_complex:
        return npy_cabs(npy_cdouble_from_double_complex(x))
    else:
        return libc.math.fabs(x)

cdef inline number_t zlog(number_t x) nogil:
    cdef np.npy_cdouble r
    if number_t is double_complex:
        r = npy_clog(npy_cdouble_from_double_complex(x))
        return double_complex_from_npy_cdouble(r)
    else:
        return libc.math.log(x)

cdef inline number_t zexp(number_t x) nogil:
    cdef np.npy_cdouble r
    if number_t is double_complex:
        r = npy_cexp(npy_cdouble_from_double_complex(x))
        return double_complex_from_npy_cdouble(r)
    else:
        return libc.math.exp(x)

