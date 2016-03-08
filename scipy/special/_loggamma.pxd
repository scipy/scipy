# An implementation of the principle branch of the logarithm of gamma
# based on the paper "Computing the Principle Branch of log-Gamma"
# (1997) by D.E.G. Hare.
#
# Author: Josh Wilson
#
# Distributed under the same license as Scipy.

import cython
from libc.math cimport M_PI, ceil, sin, log
from _complexstuff cimport nan, zisnan, zabs, zlog, zsin, zarg

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
# Euler-Mascheroni constant.
DEF EULMASCH = 0.5772156649015328606
# (EULMASCH**2 + pi**2/6)/2, the O(z) term in the Laurent series of
# the Gamma function around 0.
DEF LAURENTZ = 0.9890559953279725


@cython.cdivision(True)
cdef inline double complex loggamma(double complex z) nogil:
    """The strategy for computing loggamma is:
    - If z is in a small strip around the negative axis, use the
    reflection formula.
    - If z has negative imaginary part, conjugate it and use the
    relation loggamma(z*)* = loggamma(z)
    - If abs(z) is very small, take the log of the Laurent series for
    Gamma around the origin.
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
        return nan + 1J*nan
    elif z == 1:
        return 0

    if rz < 0 and -smallimag <= iz and iz <= smallimag:
        # Reflection formula for loggamma; see Proposition 3.1. This
        # part computes the terms log(pi/sin(pi*z)) + 2*s*k*pi*i. The
        # strategy is:
        # - If Im(z) < 0, conjugate z. This only negates the imaginary
        # part of the result.
        # - Compute pi/sin(pi*z).
        # - Compute log(abs(pi/sin(pi*z))).
        # - Compute arg(pi/sin(pi*z)).
        # - If Im(z) > 0, every point of the form m - 0.5 for even m
        # between Re(z) and 0 contributes a factor of -2*pi*i. (These
        # are the points where log(pi/sin(pi*z)) hits branch cuts.)
        # - Because of rounding errors log(pi/sin(pi*z)) might cross
        # the branch cut after z passes m - 0.5. Fix these cases.
        # - If Im(z) = 0, use Corollary 3.3.

        # Shift z to the origin before computing sin; otherwise we
        # lose accuracy around the poles.
        if iz > 0:
            logarg = M_PI/zsin(M_PI*shift(z))
        elif iz == 0:
            # If we don't treat this case seperately the imaginary
            # part will come out to be -0j, which puts us on the wrong
            # side of the branch cut in clog.
            logarg = M_PI/sin(M_PI*shift(z).real)
        else:
            logarg = M_PI/zsin(M_PI*shift(z).conjugate())
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

        if iz >= 0:
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
        if absz < 1e-4:
            # Use the Laurent expansion of Gamma. Use the bound 1e-4
            # because the singular term in the Laurent expansion is
            # O(1/z), so dropping the quadratic term should have no
            # effect to double precision.
            logterm = zlog(1.0/z - EULMASCH + LAURENTZ*z)
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


cdef inline double complex shift(double complex z) nogil:
    """Given z < 0, find shiftz so that Re(shiftz) is in [-0.5, 0.5] and
    sin(pi*shiftz) = sin(pi*z).

    """
    cdef:
        int n = <int>ceil(z.real)
        double complex shiftz

    if n % 2 == 1:
        n -= 1
    shiftz = z - n
    if shiftz.real > 0.5:
        return 1 - shiftz
    elif shiftz.real < -0.5:
        return -1 - shiftz
    else:
        return shiftz


cdef inline double find_m(double x) nogil:
    """Find the smallest m = 2*j so that x <= m - 0.5. In theory m is an
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
    """Forward recursion if s = 1 and backward recursion if s = -1. See
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
    # The Bernoulli numbers B_2k for 1 <= k <= 16.
    cdef:
        int n = 1
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
        double complex rsqz = 1/z**2
        double complex coeff = 1.0/z

    res += coeff*bernoulli2k[n-1]/(2*n*(2*n - 1))
    for n in range(2, 17):
        coeff *= rsqz
        res += coeff*bernoulli2k[n-1]/(2*n*(2*n - 1))
    return res
