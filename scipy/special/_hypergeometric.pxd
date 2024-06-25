from libc.math cimport fabs, exp, floor, isnan, M_PI, NAN, INFINITY

import cython

from . cimport sf_error

cdef extern from 'special_wrappers.h':
    double hypU_wrap(double, double, double) nogil
    double cephes_poch_wrap(double x, double m) nogil


@cython.cdivision(True)
cdef inline double hyperu(double a, double b, double x) noexcept nogil:
    if isnan(a) or isnan(b) or isnan(x):
        return NAN

    if x < 0.0:
        sf_error.error("hyperu", sf_error.DOMAIN, NULL)
        return NAN

    if x == 0.0:
        if b > 1.0:
            # DMLF 13.2.16-18
            sf_error.error("hyperu", sf_error.SINGULAR, NULL)
            return INFINITY
        else:
            # DLMF 13.2.14-15 and 13.2.19-21
            return cephes_poch_wrap(1.0 - b + a, -a)

    return hypU_wrap(a, b, x)
