
import cython

from libc.math cimport log, exp, fabs, sqrt

cdef extern from "_c99compat.h":
    int sc_isnan(double x) nogil
    int sc_isinf(double x) nogil

from ._cephes cimport ellpk
from ._complexstuff cimport pi, nan


cdef inline double _agm_iter(double a, double b) nogil:
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
cdef inline double agm(double a, double b) nogil:
    # Arithmetic-geometric mean

    # sqrthalfmax is sqrt(np.finfo(1.0).max/2)
    # invsqrthalfmax is 1/sqrthalfmax
    cdef double sqrthalfmax = 9.480751908109176e+153
    cdef double invsqrthalfmax = 1.0547686614863e-154

    cdef double e
    cdef int sgn

    if sc_isnan(a) or sc_isnan(b):
        return nan

    if (a < 0 and b > 0) or (a > 0 and b < 0):
        # a and b have opposite sign.
        return nan

    if (sc_isinf(a) or sc_isinf(b)) and (a == 0 or b == 0):
        # One value is inf and the other is 0.
        return nan

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
        return sgn*(pi/4)*(a + b)/ellpk(e)
    else:
        # At least one value is "extreme" (very big or very small).
        # Use the iteration to avoid overflow or underflow.
        return sgn*_agm_iter(a, b)
