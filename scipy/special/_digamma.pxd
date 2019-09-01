# An implementation of the digamma function for complex arguments.
#
# Author: Josh Wilson
#
# Distributed under the same license as Scipy.
#
# Sources:
# [1] "The Digital Library of Mathematical Functions", dlmf.nist.gov
#
# [2] mpmath (version 0.19), http://mpmath.org
#

import cython
from libc.math cimport ceil, fabs, M_PI
from ._complexstuff cimport number_t, nan, zlog, zabs
from ._trig cimport sinpi, cospi
from ._cephes cimport zeta, psi
from . cimport sf_error

# Use the asymptotic series for z away from the negative real axis
# with abs(z) > smallabsz.
DEF smallabsz = 16
# Use the reflection principle for z with z.real < 0 that are within
# smallimag of the negative real axis.
DEF smallimag = 6
# Relative tolerance for series
DEF tol = 2.220446092504131e-16
# All of the following were computed with mpmath
# Location of the positive root
DEF posroot = 1.4616321449683623
# Value of the positive root
DEF posrootval = -9.2412655217294275e-17
# Location of the negative root
DEF negroot = -0.504083008264455409
# Value of the negative root
DEF negrootval = 7.2897639029768949e-17


cdef inline double digamma(double z) nogil:
    """Wrap Cephes' psi to take advantage of the series expansion around
    the smallest negative zero.

    """
    if zabs(z - negroot) < 0.3:
        return zeta_series(z, negroot, negrootval)
    else:
        return psi(z)


@cython.cdivision(True)
cdef inline double complex cdigamma(double complex z) nogil:
    """
    Compute the digamma function for complex arguments. The strategy
    is:

    - Around the two zeros closest to the origin (posroot and negroot)
    use a Taylor series with precomputed zero order coefficient.
    - If close to the origin, use a recurrence relation to step away
    from the origin.
    - If close to the negative real axis, use the reflection formula
    to move to the right halfplane.
    - If |z| is large (> 16), use the asymptotic series.
    - If |z| is small, use a recurrence relation to make |z| large
    enough to use the asymptotic series.

    """
    cdef:
        int n
        double absz = zabs(z)
        double complex res = 0
        double complex init

    if z.real <= 0 and ceil(z.real) == z:
        # Poles
        sf_error.error("digamma", sf_error.SINGULAR, NULL)
        return nan + 1j*nan
    elif zabs(z - negroot) < 0.3:
        # First negative root
        return zeta_series(z, negroot, negrootval)

    if z.real < 0 and fabs(z.imag) < smallabsz:
        # Reflection formula for digamma. See
        #
        # https://dlmf.nist.gov/5.5#E4
        #
        res -= M_PI*cospi(z)/sinpi(z)
        z = 1 - z
        absz = zabs(z)

    if absz < 0.5:
        # Use one step of the recurrence relation to step away from
        # the pole.
        res -= 1/z
        z += 1
        absz = zabs(z)

    if zabs(z - posroot) < 0.5:
        res += zeta_series(z, posroot, posrootval)
    elif absz > smallabsz:
        res += asymptotic_series(z)
    elif z.real >= 0:
        n = <int>(smallabsz - absz) + 1
        init = asymptotic_series(z + n)
        res += backward_recurrence(z + n, init, n)
    else:
        # z.real < 0, absz < smallabsz, and z.imag > smallimag
        n = <int>(smallabsz - absz) - 1
        init = asymptotic_series(z - n)
        res += forward_recurrence(z - n, init, n)
    return res


@cython.cdivision(True)
cdef inline double complex forward_recurrence(double complex z,
                                              double complex psiz,
                                               int n) nogil:
    """
    Compute digamma(z + n) using digamma(z) using the recurrence
    relation

    digamma(z + 1) = digamma(z) + 1/z.

    See https://dlmf.nist.gov/5.5#E2

    """
    cdef:
        int k
        double complex res = psiz

    for k in range(n):
        res += 1/(z + k)
    return res


@cython.cdivision(True)
cdef inline double complex backward_recurrence(double complex z,
                                               double complex psiz,
                                               int n) nogil:
    """
    Compute digamma(z - n) using digamma(z) and a recurrence
    relation.

    """
    cdef:
        int k
        double complex res = psiz

    for k in range(1, n + 1):
        res -= 1/(z - k)
    return res


@cython.cdivision(True)
cdef inline double complex asymptotic_series(double complex z) nogil:
    """
    Evaluate digamma using an asymptotic series. See

    https://dlmf.nist.gov/5.11#E2

    """
    cdef:
        int k = 1
        # The Bernoulli numbers B_2k for 1 <= k <= 16.
        double *bernoulli2k = [
            0.166666666666666667, -0.0333333333333333333,
            0.0238095238095238095, -0.0333333333333333333,
            0.0757575757575757576, -0.253113553113553114,
            1.16666666666666667, -7.09215686274509804,
            54.9711779448621554, -529.124242424242424,
            6192.12318840579710, -86580.2531135531136,
            1425517.16666666667, -27298231.0678160920,
            601580873.900642368, -15116315767.0921569]
        double complex rzz = 1/z/z
        double complex zfac = 1
        double complex term
        double complex res

    res = zlog(z) - 0.5/z
    for k in range(1, 17):
        zfac *= rzz
        term = -bernoulli2k[k-1]*zfac/(2*k)
        res += term
        if zabs(term) < tol*zabs(res):
            break
    return res


cdef inline number_t zeta_series(number_t z, double root, double rootval) nogil:
    """
    The coefficients of the Taylor series for digamma at any point can
    be expressed in terms of the Hurwitz zeta function. If we
    precompute the floating point number closest to a zero and the 0th
    order Taylor coefficient at that point then we can compute higher
    order coefficients without loss of accuracy using zeta (the zeros
    are simple) and maintain high-order accuracy around the zeros.

    """
    cdef:
        int n
        number_t res = rootval
        number_t coeff = -1
        number_t term

    z = z - root
    for n in range(1, 100):
        coeff *= -z
        term = coeff*zeta(n + 1, root)
        res += term
        if zabs(term) < tol*zabs(res):
            break
    return res
