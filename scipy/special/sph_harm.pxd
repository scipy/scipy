from . cimport sf_error

cdef extern from "specfun_wrappers.h":
    double pmv_wrap(double, double, double) nogil

cdef extern from "c_misc/misc.h":
    double poch(double x, double m) nogil

from ._complexstuff cimport *
from libc.math cimport cos, sin, sqrt, fabs, sqrt, exp, M_PI
from libc.stdlib cimport abs

cdef inline double complex sph_harmonic(int m, int n, double theta, double phi) nogil:
    cdef double x, prefactor, k
    cdef double complex val, l
    cdef int mp

    if m == n and m > 65:
        # Fall back to a sketchy implementation that returns
        # a somewhat reasonable result, as the correct implementation fails
        k = (m*phi)
        l = cos(k) + 1j*sin(k)
        val = sqrt(1/M_PI)*0.5524*m**0.2428*l*(sin(theta))**m

        if val.real > 1E-7 or val.imag > 1E-7:
            return val

    x = cos(phi)
    if abs(m) > n :
        sf_error.error("sph_harm", sf_error.ARG, "m should not be greater than n")
        return nan
    if n < 0:
        sf_error.error("sph_harm", sf_error.ARG, "n should not be negative")
        return nan
    if m < 0:
        mp = -m
        prefactor = (-1)**mp * poch(n + mp + 1, -2 * mp)
    else:
        mp = m
    val = pmv_wrap(mp, n, x)
    if  m < 0:
        val *= prefactor
    val *= sqrt((2*n + 1) / 4.0 / pi)
    val *= sqrt(poch(n + m + 1, -2 * m))
    val *= zexp(1j * m * theta)
    return val
