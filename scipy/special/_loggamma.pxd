# An implementation of the principal branch of the logarithm of
# Gamma. Also contains implementations of Gamma and 1/Gamma which are
# easily computed from log-Gamma.
#
# Author: Josh Wilson
#
# Distributed under the same license as Scipy.
#
# References
# ----------
# [1] Hare, "Computing the Principal Branch of log-Gamma",
#     Journal of Algorithms, 1997.
#
# [2] Julia,
#     https://github.com/JuliaLang/julia/blob/master/base/special/gamma.jl
#
cimport cython
from . cimport sf_error
from libc.math cimport M_PI, NAN, floor, fabs, copysign, signbit

from ._complexstuff cimport (
    zisnan, zabs, zlog, zlog1, zexp, zpack
)
from ._trig cimport sinpi
from ._evalpoly cimport cevalpoly

cdef extern from "cephes.h":
    double lgam(double x) nogil

DEF TWOPI = 6.2831853071795864769252842 # 2*pi
DEF LOGPI = 1.1447298858494001741434262 # log(pi)
DEF HLOG2PI = 0.918938533204672742 # log(2*pi)/2
DEF SMALLX = 7
DEF SMALLY = 7
DEF TAYLOR_RADIUS = 0.2


cdef inline double loggamma_real(double x) noexcept nogil:
    if x < 0.0:
        return NAN

    return lgam(x)


@cython.cdivision(True)
cdef inline double complex loggamma(double complex z) noexcept nogil:
    """Compute the principal branch of log-Gamma."""
    cdef double tmp

    if zisnan(z):
        return zpack(NAN, NAN)
    elif z.real <= 0 and z == floor(z.real):
        sf_error.error("loggamma", sf_error.SINGULAR, NULL)
        return zpack(NAN, NAN)
    elif z.real > SMALLX or fabs(z.imag) > SMALLY:
        return loggamma_stirling(z)
    elif zabs(z - 1) <= TAYLOR_RADIUS:
        return loggamma_taylor(z)
    elif zabs(z - 2) <= TAYLOR_RADIUS:
        # Recurrence relation and the Taylor series around 1
        return zlog1(z - 1) + loggamma_taylor(z - 1)
    elif z.real < 0.1:
        # Reflection formula; see Proposition 3.1 in [1]
        tmp = copysign(TWOPI, z.imag)*floor(0.5*z.real + 0.25)
        return zpack(LOGPI, tmp) - zlog(sinpi(z)) - loggamma(1 - z)
    elif signbit(z.imag) == 0:
        # z.imag >= 0 but is not -0.0
        return loggamma_recurrence(z)
    else:
        return loggamma_recurrence(z.conjugate()).conjugate()


@cython.cdivision(True)
cdef inline double complex loggamma_recurrence(double complex z) noexcept nogil:
    """Backward recurrence relation.

    See Proposition 2.2 in [1] and the Julia implementation [2].

    """
    cdef:
        int signflips = 0
        int sb = 0
        int nsb
        double complex shiftprod = z

    z.real += 1
    while z.real <= SMALLX:
        shiftprod *= z
        nsb = signbit(shiftprod.imag)
        signflips += 1 if nsb != 0 and sb == 0 else 0
        sb = nsb
        z.real += 1
    return loggamma_stirling(z) - zlog(shiftprod) - signflips*TWOPI*1J


@cython.cdivision(True)
cdef inline double complex loggamma_stirling(double complex z) noexcept nogil:
    """Stirling series for log-Gamma.

    The coefficients are B[2*n]/(2*n*(2*n - 1)) where B[2*n] is the
    (2*n)th Bernoulli number. See (1.1) in [1].

    """
    cdef:
        double *coeffs = [
            -2.955065359477124183e-2, 6.4102564102564102564e-3,
            -1.9175269175269175269e-3, 8.4175084175084175084e-4,
            -5.952380952380952381e-4, 7.9365079365079365079e-4,
            -2.7777777777777777778e-3, 8.3333333333333333333e-2
        ]
        double complex rz = 1.0/z
        double complex rzz = rz/z

    return (z - 0.5)*zlog(z) - z + HLOG2PI + rz*cevalpoly(coeffs, 7, rzz)


@cython.cdivision(True)
cdef inline double complex loggamma_taylor(double complex z) noexcept nogil:
    """Taylor series for log-Gamma around z = 1.

    It is

    loggamma(z + 1) = -gamma*z + zeta(2)*z**2/2 - zeta(3)*z**3/3 ...

    where gamma is the Euler-Mascheroni constant.

    """
    cdef:
        double *coeffs = [
            -4.3478266053040259361e-2, 4.5454556293204669442e-2,
            -4.7619070330142227991e-2, 5.000004769810169364e-2,
            -5.2631679379616660734e-2, 5.5555767627403611102e-2,
            -5.8823978658684582339e-2, 6.2500955141213040742e-2,
            -6.6668705882420468033e-2, 7.1432946295361336059e-2,
            -7.6932516411352191473e-2, 8.3353840546109004025e-2,
            -9.0954017145829042233e-2, 1.0009945751278180853e-1,
            -1.1133426586956469049e-1, 1.2550966952474304242e-1,
            -1.4404989676884611812e-1, 1.6955717699740818995e-1,
            -2.0738555102867398527e-1, 2.7058080842778454788e-1,
            -4.0068563438653142847e-1, 8.2246703342411321824e-1,
            -5.7721566490153286061e-1
        ]

    z = z - 1
    return z*cevalpoly(coeffs, 22, z)


cdef inline double complex cgamma(double complex z) noexcept nogil:
    """Compute Gamma(z) using loggamma."""
    if z.real <= 0 and z == floor(z.real):
        # Poles
        sf_error.error("gamma", sf_error.SINGULAR, NULL)
        return NAN + 1J*NAN
    return zexp(loggamma(z))


cdef inline double complex crgamma(double complex z) noexcept nogil:
    """Compute 1/Gamma(z) using loggamma."""
    if z.real <= 0 and z == floor(z.real):
        # Zeros at 0, -1, -2, ...
        return 0
    return zexp(-loggamma(z))
