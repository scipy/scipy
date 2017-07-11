"""
Core functions for the Frechet distribution ufuncs.
"""

from libc.math cimport log, log1p, expm1, exp, pow
from ._complexstuff cimport inf


cdef inline double _frechet_logpdf(double x, double alpha) nogil:
    if x <= 0:
        return -inf
    else:
        return log(alpha) + (-alpha-1)*log(x) - pow(x, -alpha)


cdef inline double _frechet_pdf(double x, double alpha) nogil:  
    if x <= 0:
        return 0.0
    else:
        return exp(_frechet_logpdf(x, alpha))


cdef inline double _frechet_cdf(double x, double alpha) nogil:
    if  x <= 0:
        return 0.0
    else:
        return exp(-pow(x, -alpha))


cdef inline double _frechet_ppf(double p, double alpha) nogil:
    if p == 1.0:
        return inf
    elif p == 0.0:
        return 0.0
    else:
        return pow(-1/log(p), 1/alpha)


cdef inline double _frechet_sf(double x, double alpha) nogil:
    if x <= 0:
        return 1.0
    else:
        return -expm1(-pow(x, -alpha))


cdef inline double _frechet_isf(double p, double alpha) nogil:
    if p == 1.0:
        return 0.0
    elif p == 0.0:
        return inf
    else:
        return pow(-1/log1p(-p), 1/alpha)

