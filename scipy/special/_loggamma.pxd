# An implementation of the principal branch of the logarithm of gamma
# based on the paper "Computing the Principal Branch of log-Gamma"
# (1997) by D.E.G. Hare. Also contains implementations of gamma and
# 1/gamma which are easily computed from loggamma.
#
# Author: Josh Wilson
#
# Distributed under the same license as Scipy.

import cython
cimport sf_error
from libc.math cimport M_PI, ceil, sin, log
from _complexstuff cimport (
    nan, zisnan, zabs, zlog, zlog1, zsin, zarg, zexp, zlog, zdiv
)
from _trig cimport sinpi

cdef extern from "numpy/npy_math.h":
    int npy_signbit(double x) nogil
    double NPY_EULER

cdef extern from "cephes.h":
    double zeta(double x, double q) nogil

# log(2*pi)/2
DEF HLOG2PI = 0.918938533204672742
# Use the recurrence relation to make |z| bigger than smallz before
# using the asymptotic series. See (4.3) for information about picking
# it.
DEF smallz = 16
# Can use the asymptotic series in the left halfplane if the imaginary
# part of z is above this value. The number 16 comes from the number
# of digits of desired precision; see (4.4).
DEF smallimag = 0.37*16
# Relative tolerance for series expansions
DEF tol = 2.2204460492503131e-16


@cython.cdivision(True)
cdef inline double complex loggamma(double complex z) nogil:
    """
    The strategy for computing loggamma is:
    - If z is in a small strip around the negative axis, use the
    reflection formula.
    - If z has negative imaginary part, conjugate it and use the
    relation loggamma(z*)* = loggamma(z)
    - If z is close to the zero at 1, use a Taylor series.
    - If z is close to 0 or 2, use a recurence relation and the Taylor
    series at 1.
    - If abs(z) is large, use an asymptotic series.
    - Else use a recurrence relation to make z larger and then use
    the asymptotic series.

    """
    cdef:
        int conjugated = 0
        int reflected = 0
        int n
        double rz = z.real
        double iz = z.imag
        double absz = zabs(z)
        double m
        double complex res = 0
        double complex init, logarg, argterm, logterm

    if zisnan(z):
        return z
    elif rz <= 0 and z == ceil(rz):
        # Poles
        sf_error.error("loggamma", sf_error.SINGULAR, NULL)
        return nan + 1J*nan

    if rz < 0 and -smallimag <= iz and iz <= smallimag:
        # Reflection formula for loggamma; see Proposition 3.1. This
        # part computes the terms log(pi/sin(pi*z)) + 2*s*k*pi*i. The
        # strategy is:
        # - If Im(z) < 0, conjugate z. This only negates the imaginary
        # part of the result.
        # - Compute pi/sin(pi*z).
        # - Compute real and imaginary parts of log(pi/sin(pi*z)).
        # - If Im(z) > 0, every point of the form m - 0.5 for even m
        # between Re(z) and 0 contributes a factor of -2*pi*i. (These
        # are the points where log(pi/sin(pi*z)) hits branch cuts.)
        # - Because of rounding errors log(pi/sin(pi*z)) might cross
        # the branch cut after z passes m - 0.5. Fix these cases.
        # - If Im(z) = 0, use Corollary 3.3.
        if iz > 0:
            logarg = M_PI/sinpi(z)
        elif iz == 0:
            # If we don't treat this case seperately the imaginary
            # part will come out to be -0j, which puts us on the wrong
            # side of the branch cut in clog.
            logarg = M_PI/sinpi(z.real)
        else:
            logarg = M_PI/sinpi(z.conjugate())
        res += log(zabs(logarg))
        argterm = zarg(logarg)

        if iz == 0:
            argterm += 2*M_PI*ceil(rz/2 - 1)
        elif rz <= -0.5:
            m = find_m(rz)
            argterm += (m - 2)*M_PI
            if rz > m - 1.5 and logarg.real < 0 and logarg.imag < 0:
                # rz crossed the branch cut but due to rounding errors
                # log(pi/sin(pi*z)) didn't.
                argterm += 2*M_PI

        if npy_signbit(iz) == 0:
            res += 1j*argterm
        else:
            res -= 1j*argterm
        z = 1 - z
        rz, iz, absz = z.real, z.imag, zabs(z)
        reflected = 1
    if iz < 0:
        # loggamma(z) = loggamma(z*)*
        z = z.conjugate()
        rz, iz, absz = z.real, z.imag, zabs(z)
        conjugated = 1

    if rz >= 0:
        if zabs(z - 1) <= 0.5:
            logterm = taylor(z)
        elif zabs(z - 2) < 0.5:
            # Use the recurrence relation and Taylor series around 1.
            logterm = zlog1(z - 1) + taylor(z - 1)
        elif absz < 0.5:
            # Use the recurrence relation and Taylor series around 1.
            logterm = -zlog(z) + taylor(z + 1)
        elif absz >= smallz:
            logterm = asymptotic_series(z)
        else:
            n = <int>ceil(smallz - rz)
            init = asymptotic_series(z + n)
            logterm = recurrence(z + n, init, n, -1)
    else:
        # Case where creal(z) < 0 and cimag(z) > smallimag.
        if absz >= smallz:
            logterm = asymptotic_series(z)
        else:
            n = <int>ceil(smallz + rz)
            init = asymptotic_series(z - n)
            logterm = recurrence(z - n, init, n, 1)

    if conjugated:
        logterm = logterm.conjugate()
    if reflected:
        res -= logterm
    else:
        res = logterm
    return res


cdef inline double find_m(double x) nogil:
    """
    Find the smallest m = 2*j so that x <= m - 0.5. In theory m is an
    integer, make it a double to avoid overflow.

    """
    cdef:
        double m = ceil(x)
        double hm = m/2

    if hm == ceil(hm):
        if m - x >= 0.5:
            return m
        else:
            return m + 2
    else:
        return m + 1


cdef inline int imag_sgncmp(double complex z1, double complex z2) nogil:
    return 1 if z1.imag >= 0 and z2.imag < 0 else 0


@cython.cdivision(True)
cdef inline double complex recurrence(double complex z, double complex
                                      init, int n, int s) nogil:
    """
    Forward recursion if s = 1 and backward recursion if s = -1. See
    Proposition 2.2.

    """
    cdef:
        int k = 0
        int i
        double complex po, pn

    if s == 1:
        po = z
        # Make sure pn is set even if we don't go through the loop.
        pn = po
        for i in range(1, n):
            pn = po*(z + i)
            k += imag_sgncmp(po, pn)
            po = pn
    else:
        po = z - 1
        pn = po
        for i in range(2, n + 1):
            pn = po*(z - i)
            k += imag_sgncmp(po, pn)
            po = pn
    return init + s*(zlog(pn) + 2*k*M_PI*1J)


@cython.cdivision(True)
cdef inline double complex asymptotic_series(double complex z) nogil:
    cdef:
        int n = 1
        # The Bernoulli numbers B_2k for 1 <= k <= 16.
        double *bernoulli2k = [0.166666666666666667,
                               -0.0333333333333333333,
                               0.0238095238095238095,
                               -0.0333333333333333333,
                               0.0757575757575757576,
                               -0.253113553113553114,
                               1.16666666666666667,
                               -7.09215686274509804,
                               54.9711779448621554,
                               -529.124242424242424,
                               6192.12318840579710,
                               -86580.2531135531136,
                               1425517.16666666667,
                               -27298231.0678160920,
                               601580873.900642368,
                               -15116315767.0921569]
        double complex res = (z - 0.5)*zlog(z) - z + HLOG2PI
        double complex coeff = zdiv(1.0, z)
        double complex rsqz = zdiv(coeff, z)
        double complex term

    res += coeff*bernoulli2k[n-1]/(2*n*(2*n - 1))
    for n in range(2, 17):
        coeff *= rsqz
        term = coeff*bernoulli2k[n-1]/(2*n*(2*n - 1))
        res += term
        if zabs(term) <= tol*zabs(res):
            break
    return res


@cython.cdivision(True)
cdef inline double complex taylor(double complex z) nogil:
    """
    Taylor series around z = 1. Derived from the Taylor series for the
    digamma function at 1:

    psi(z + 1) = -gamma - (-zeta(2)*z + zeta(3)*z**2 - ...)

    by integrating. Here gamma is the Euler-Mascheroni constant and
    zeta is the Riemann zeta function. Note that we can approximate
    zeta(n) by 1/(n - 1), so if we use a radius of 1/2 we should need
    at most 41 terms.

    """
    cdef:
        int n
        double complex zfac, coeff, res

    z = z - 1
    if z == 0:
        return 0
    res = -NPY_EULER*z
    zfac = -z
    for n in xrange(2, 42):
        zfac *= -z
        coeff = zeta(n, 1)*zfac/n
        res += coeff
        if zabs(coeff/res) < tol:
            break
    return res


cdef inline double complex cgamma(double complex z) nogil:
    """Compute Gamma(z) using loggamma."""
    if z.real <= 0 and z == ceil(z.real):
        # Poles
        sf_error.error("gamma", sf_error.SINGULAR, NULL)
        return nan + 1J*nan
    return zexp(loggamma(z))


cdef inline double complex crgamma(double complex z) nogil:
    """Compute 1/Gamma(z) using loggamma."""
    if z.real <= 0 and z == ceil(z.real):
        # Zeros at 0, -1, -2, ...
        return 0
    return zexp(-loggamma(z))
