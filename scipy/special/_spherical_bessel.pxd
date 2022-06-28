# -*-cython-*-
#
# Implementation of spherical Bessel functions and modified spherical Bessel
# functions of the first and second kinds, as well as their derivatives.
#
# Author: Tadeusz Pudlik
#
# Distributed under the same license as SciPy.
#
# I attempt to correctly handle the edge cases (0 and infinity), but this is
# tricky: the values of the functions often depend on the direction in which
# the limit is taken. At zero, I follow the convention of numpy (1.9.2),
# which treats zero differently depending on its type:
#
#   >>> np.cos(0)/0
#   inf
#   >>> np.cos(0+0j)/(0+0j)
#   inf + nan*j
#
# So, real zero is assumed to be "positive zero", while complex zero has an
# unspecified sign and produces nans.  Similarly, complex infinity is taken to
# represent the "point at infinity", an ambiguity which for some functions
# makes `nan` the correct return value.

import cython
from libc.math cimport cos, sin, sqrt, M_PI_2

from numpy cimport npy_cdouble
from ._complexstuff cimport *

from . cimport sf_error

cdef extern from "amos_wrappers.h":
    npy_cdouble cbesi_wrap( double v, npy_cdouble z) nogil
    npy_cdouble cbesj_wrap(double v, npy_cdouble z) nogil
    double cbesj_wrap_real(double v, double x) nogil
    npy_cdouble cbesk_wrap(double v, npy_cdouble z) nogil
    double cbesk_wrap_real(double v, double x) nogil
    npy_cdouble cbesy_wrap(double v, npy_cdouble z) nogil
    double cbesy_wrap_real(double v, double x) nogil

from ._cephes cimport iv

# Fused type wrappers

cdef inline number_t cbesj(double v, number_t z) nogil:
    cdef npy_cdouble r
    if number_t is double:
        return cbesj_wrap_real(v, z)
    else:
        r = cbesj_wrap(v, npy_cdouble_from_double_complex(z))
        return double_complex_from_npy_cdouble(r)

cdef inline number_t cbesy(double v, number_t z) nogil:
    cdef npy_cdouble r
    if number_t is double:
        return cbesy_wrap_real(v, z)
    else:
        r = cbesy_wrap(v, npy_cdouble_from_double_complex(z))
        return double_complex_from_npy_cdouble(r)

cdef inline number_t cbesk(double v, number_t z) nogil:
    cdef npy_cdouble r
    if number_t is double:
        return cbesk_wrap_real(v, z)
    else:
        r = cbesk_wrap(v, npy_cdouble_from_double_complex(z))
        return double_complex_from_npy_cdouble(r)


# Spherical Bessel functions

@cython.cdivision(True)
cdef inline double spherical_jn_real(long n, double x) nogil:
    cdef double s0, s1, sn
    cdef int idx

    if npy_isnan(x):
        return x
    if n < 0:
        sf_error.error("spherical_jn", sf_error.DOMAIN, NULL)
        return nan
    if x == inf or x == -inf:
        return 0
    if x == 0:
        if n == 0:
            return 1
        else:
            return 0

    if n > 0 and n >= x:
        return sqrt(M_PI_2/x)*cbesj(n + 0.5, x)

    s0 = sin(x)/x
    if n == 0:
        return s0
    s1 = (s0 - cos(x))/x
    if n == 1:
        return s1

    for idx in range(n - 1):
        sn = (2*idx + 3)*s1/x - s0
        s0 = s1
        s1 = sn
        if npy_isinf(sn):
            # Overflow occurred already: terminate recurrence.
            return sn

    return sn


@cython.cdivision(True)
cdef inline double complex spherical_jn_complex(long n, double complex z) nogil:
    cdef double complex out
    if zisnan(z):
        return z
    if n < 0:
        sf_error.error("spherical_jn", sf_error.DOMAIN, NULL)
        return nan
    if z.real == inf or z.real == -inf:
        # https://dlmf.nist.gov/10.52.E3
        if z.imag == 0:
            return 0
        else:
            return (1+1j)*inf
    if z.real == 0 and z.imag == 0:
        if n == 0:
            return 1
        else:
            return 0

    out = zsqrt(M_PI_2/z)*cbesj(n + 0.5, z)

    if z.imag == 0:
        # Small imaginary part is spurious
        return out.real
    else:
        return out


@cython.cdivision(True)
cdef inline double spherical_yn_real(long n, double x) nogil:
    cdef double s0, s1, sn
    cdef int idx

    if npy_isnan(x):
        return x
    if n < 0:
        sf_error.error("spherical_yn", sf_error.DOMAIN, NULL)
        return nan
    if x < 0:
        return (-1)**(n+1)*spherical_yn_real(n, -x)
    if x == inf or x == -inf:
        return 0
    if x == 0:
        return -inf

    s0 = -cos(x)/x
    if n == 0:
        return s0
    s1 = (s0 - sin(x))/x
    if n == 1:
        return s1

    for idx in range(n - 1):
        sn = (2*idx + 3)*s1/x - s0
        s0 = s1
        s1 = sn
        if npy_isinf(sn):
            # Overflow occurred already: terminate recurrence.
            return sn

    return sn


@cython.cdivision(True)
cdef inline double complex spherical_yn_complex(long n, double complex z) nogil:

    if zisnan(z):
        return z
    if n < 0:
        sf_error.error("spherical_yn", sf_error.DOMAIN, NULL)
        return nan
    if z.real == 0 and z.imag == 0:
        # https://dlmf.nist.gov/10.52.E2
        return nan
    if z.real == inf or z.real == -inf:
        # https://dlmf.nist.gov/10.52.E3
        if z.imag == 0:
            return 0
        else:
            return (1+1j)*inf

    return zsqrt(M_PI_2/z)*cbesy(n + 0.5, z)


@cython.cdivision(True)
cdef inline double spherical_in_real(long n, double z) nogil:

    if npy_isnan(z):
        return z
    if n < 0:
        sf_error.error("spherical_in", sf_error.DOMAIN, NULL)
        return nan
    if z == 0:
        # https://dlmf.nist.gov/10.52.E1
        if n == 0:
            return 1
        else:
            return 0
    if npy_isinf(z):
        # https://dlmf.nist.gov/10.49.E8
        if z == -inf:
            return (-1)**n*inf
        else:
            return inf

    return sqrt(M_PI_2/z)*iv(n + 0.5, z)


@cython.cdivision(True)
cdef inline double complex spherical_in_complex(long n, double complex z) nogil:
    cdef npy_cdouble s

    if zisnan(z):
        return z
    if n < 0:
        sf_error.error("spherical_in", sf_error.DOMAIN, NULL)
        return nan
    if zabs(z) == 0:
        # https://dlmf.nist.gov/10.52.E1
        if n == 0:
            return 1
        else:
            return 0
    if zisinf(z):
        # https://dlmf.nist.gov/10.52.E5
        if z.imag == 0:
            if z.real == -inf:
                return (-1)**n*inf
            else:
                return inf
        else:
            return nan

    s = cbesi_wrap(n + 0.5, npy_cdouble_from_double_complex(z))
    return zsqrt(M_PI_2/z)*double_complex_from_npy_cdouble(s)


@cython.cdivision(True)
cdef inline double spherical_kn_real(long n, double z) nogil:

    if npy_isnan(z):
        return z
    if n < 0:
        sf_error.error("spherical_kn", sf_error.DOMAIN, NULL)
        return nan
    if z == 0:
        return inf
    if npy_isinf(z):
        # https://dlmf.nist.gov/10.52.E6
        if z == inf:
            return 0
        else:
            return -inf

    return sqrt(M_PI_2/z)*cbesk(n + 0.5, z)


@cython.cdivision(True)
cdef inline double complex spherical_kn_complex(long n, double complex z) nogil:

    if zisnan(z):
        return z
    if n < 0:
        sf_error.error("spherical_kn", sf_error.DOMAIN, NULL)
        return nan
    if zabs(z) == 0:
        return nan
    if zisinf(z):
        # https://dlmf.nist.gov/10.52.E6
        if z.imag == 0:
            if z.real == inf:
                return 0
            else:
                return -inf
        else:
            return nan

    return zsqrt(M_PI_2/z)*cbesk(n + 0.5, z)


# Derivatives

@cython.cdivision(True)
cdef inline double spherical_jn_d_real(long n, double x) nogil:
    if n == 0:
        return -spherical_jn_real(1, x)
    else:
        if x == 0:
            # DLMF 10.51.2 doesn't work, so use 10.51.1 to get the
            # exact value
            if n == 1:
                return 1.0/3
            else:
                return 0
        # DLMF 10.51.2
        return (spherical_jn_real(n - 1, x) -
                (n + 1)*spherical_jn_real(n, x)/x)


@cython.cdivision(True)
cdef inline double complex spherical_jn_d_complex(long n, double complex x) nogil:
    if n == 0:
        return -spherical_jn_complex(1, x)
    else:
        return (spherical_jn_complex(n - 1, x) -
                (n + 1)*spherical_jn_complex(n, x)/x)


@cython.cdivision(True)
cdef inline double spherical_yn_d_real(long n, double x) nogil:
    if n == 0:
        return -spherical_yn_real(1, x)
    else:
        return (spherical_yn_real(n - 1, x) -
                (n + 1)*spherical_yn_real(n, x)/x)


@cython.cdivision(True)
cdef inline double complex spherical_yn_d_complex(long n, double complex x) nogil:
    if n == 0:
        return -spherical_yn_complex(1, x)
    else:
        return (spherical_yn_complex(n - 1, x) -
                (n + 1)*spherical_yn_complex(n, x)/x)


@cython.cdivision(True)
cdef inline double spherical_in_d_real(long n, double x) nogil:
    if n == 0:
        return spherical_in_real(1, x)
    else:
        if x == 0:
            return 0
        return (spherical_in_real(n - 1, x) -
                (n + 1)*spherical_in_real(n, x)/x)


@cython.cdivision(True)
cdef inline double complex spherical_in_d_complex(long n, double complex x) nogil:
    if n == 0:
        return spherical_in_complex(1, x)
    else:
        if x == 0:
            return 0
        return (spherical_in_complex(n - 1, x) -
                (n + 1)*spherical_in_complex(n, x)/x)


@cython.cdivision(True)
cdef inline double spherical_kn_d_real(long n, double x) nogil:
    if n == 0:
        return -spherical_kn_real(1, x)
    else:
        return (-spherical_kn_real(n - 1, x) -
                (n + 1)*spherical_kn_real(n, x)/x)


@cython.cdivision(True)
cdef inline double complex spherical_kn_d_complex(long n, double complex x) nogil:
    if n == 0:
        return -spherical_kn_complex(1, x)
    else:
        return (-spherical_kn_complex(n - 1, x) -
                (n + 1)*spherical_kn_complex(n, x)/x)
