from libc.math cimport fabs, exp, floor, isnan, M_PI, NAN, INFINITY

import cython

from . cimport sf_error

cdef extern from 'xsf_wrappers.h':
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
    if b == 1 and x < 1 and -0.25 < a < 0.3:
        # DLMF 13.3.7. Fixes gh-15650
        return (x + 1 + 2*a)*hypU_wrap(a + 1, 1, x) - (a + 1)**2 * hypU_wrap(a + 2, 1, x)
    return hypU_wrap(a, b, x)
