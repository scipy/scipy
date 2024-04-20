# cython: cpow=True

import cython

from libc.math cimport log, exp, fabs, sqrt, isnan, isinf, NAN, M_PI


cdef extern from "special_wrappers.h" nogil:
    double cephes_ellpk_wrap(double x)


cdef inline double _agm_iter(double a, double b) noexcept nogil:
    # Arithmetic-geometric mean, iterative implementation
    # a and b must be positive (not zero, not nan).

    cdef double amean, gmean

    cdef int count = 20
    amean = 0.5*a + 0.5*b
    while (count > 0) and (amean != a and amean != b):
        gmean = sqrt(a)*sqrt(b)
        a = amean
        b = gmean
        amean = 0.5*a + 0.5*b
        count -= 1
    return amean


@cython.cdivision(True)
cdef inline double agm(double a, double b) noexcept nogil:
    # Arithmetic-geometric mean

    # sqrthalfmax is sqrt(np.finfo(1.0).max/2)
    # invsqrthalfmax is 1/sqrthalfmax
    cdef double sqrthalfmax = 9.480751908109176e+153
    cdef double invsqrthalfmax = 1.0547686614863e-154

    cdef double e
    cdef int sgn

    if isnan(a) or isnan(b):
        return NAN

    if (a < 0 and b > 0) or (a > 0 and b < 0):
        # a and b have opposite sign.
        return NAN

    if (isinf(a) or isinf(b)) and (a == 0 or b == 0):
        # One value is inf and the other is 0.
        return NAN

    if a == 0 or b == 0:
        # At least one of the arguments is 0.
        return 0.0

    if a == b:
        return a

    sgn = 1
    if a < 0:
        sgn = -1
        a = -a
        b = -b

    # At this point, a and b are both positive and not nan.

    if (invsqrthalfmax < a < sqrthalfmax) and (invsqrthalfmax < b < sqrthalfmax):
        e = 4*a*b/(a + b)**2
        return sgn*(M_PI/4)*(a + b)/cephes_ellpk_wrap(e)
    else:
        # At least one value is "extreme" (very big or very small).
        # Use the iteration to avoid overflow or underflow.
        return sgn*_agm_iter(a, b)
