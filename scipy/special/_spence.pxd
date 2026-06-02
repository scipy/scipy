# cython: cpow=True

# Implement Spence's function, a.k.a. the dilogarithm, for complex
# arguments. Note that our definition differs from that in the sources
# by the mapping z -> 1 - z.
#
# Sources
# [1] Zagier, "The Dilogarithm Function"
# [2] functions.wolfram.com
# [3] Ginsberg, Zaborowski, "The Dilogarithm Function of a Real Argument"
#
# Author: Josh Wilson
#
# Released under the same license as Scipy.

import cython
from ._complexstuff cimport zlog1, zabs

# Relative tolerance for the series
DEF TOL = 2.220446092504131e-16
DEF PISQ_6 = 1.6449340668482264365


@cython.cdivision(True)
cdef inline double complex cspence(double complex z) noexcept nogil:
    """
    Compute Spence's function for complex arguments. The strategy is:
    - If z is close to 0, use a series centered at 0.
    - If z is far away from 1, use the reflection formula

    spence(z) = -spence(z/(z - 1)) - pi**2/6 - ln(z - 1)**2/2

    to move close to 1. See [1].
    - If z is close to 1, use a series centered at 1.

    """
    if zabs(z) < 0.5:
        # This step isn't necessary, but this series converges faster.
        return cspence_series0(z)
    elif zabs(1 - z) > 1:
        return -cspence_series1(z/(z - 1)) - PISQ_6 - 0.5*zlog1(z - 1)**2
    else:
        return cspence_series1(z)


@cython.cdivision(True)
cdef inline double complex cspence_series0(double complex z) noexcept nogil:
    """
    A series centered at z = 0; see

    http://functions.wolfram.com/10.07.06.0005.02

    """
    cdef:
        int n
        double complex zfac = 1
        double complex sum1 = 0
        double complex sum2 = 0
        double complex term1, term2

    if z == 0:
        return PISQ_6

    for n in range(1, 500):
        zfac *= z
        term1 = zfac/n**2
        sum1 += term1
        term2 = zfac/n
        sum2 += term2
        if zabs(term1) <= TOL*zabs(sum1) and zabs(term2) <= TOL*zabs(sum2):
            break
    return PISQ_6 - sum1 + zlog1(z)*sum2


@cython.cdivision(True)
cdef inline double complex cspence_series1(double complex z) noexcept nogil:
    """
    A series centered at z = 1 which enjoys faster convergence than
    the Taylor series. See [3]. The number of terms used comes from
    bounding the absolute tolerance at the edge of the radius of
    convergence where the sum is O(1).

    """
    cdef:
        int n
        double complex zfac = 1
        double complex res = 0
        double complex term, zz

    if z == 1:
        return 0
    z = 1 - z
    zz = z**2
    for n in range(1, 500):
        zfac *= z
        # Do the divisions one at a time to guard against overflow
        term = ((zfac/n**2)/(n + 1)**2)/(n + 2)**2
        res += term
        if zabs(term) <= TOL*zabs(res):
            break
    res *= 4*zz
    res += 4*z + 5.75*zz + 3*(1 - zz)*zlog1(1 - z)
    res /= 1 + 4*z + zz
    return res
