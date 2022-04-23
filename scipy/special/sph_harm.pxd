from . cimport sf_error
from ._cephes cimport poch

cdef extern from "specfun_wrappers.h":
    double pmv_wrap(double, double, double) nogil

from ._complexstuff cimport *
from libc.math cimport cos, sqrt, fabs
from libc.stdlib cimport abs

cdef inline double complex sph_harmonic(int m, int n, double theta, double phi) nogil:
    cdef double x, prefactor
    cdef double complex val
    cdef int mp
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
