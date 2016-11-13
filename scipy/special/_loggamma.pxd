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
cimport sf_error
from libc.math cimport M_PI, floor, fabs
from _complexstuff cimport (
    nan, zisnan, zabs, zlog, zlog1, zsin, zexp, zdiv, zpack
)
from _trig cimport sinpi

cdef extern from "numpy/npy_math.h":
    double npy_copysign(double x, double y) nogil
    int npy_signbit(double x) nogil
    double NPY_EULER

cdef extern from "cephes.h":
    double zeta(double x, double q) nogil

DEF TWOPI = 6.2831853071795864769252842 # 2*pi
DEF LOGPI = 1.1447298858494001741434262 # log(pi)
DEF HLOG2PI = 0.918938533204672742 # log(2*pi)/2
DEF SMALLX = 7
DEF SMALLY = 7
DEF TOL = 2.2204460492503131e-16


@cython.cdivision(True)
cdef inline double complex loggamma(double complex z) nogil:
    """Compute the principal branch of log-Gamma."""
    cdef double tmp

    if zisnan(z):
        return zpack(nan, nan)
    elif z.real <= 0 and z == floor(z.real):
        sf_error.error("loggamma", sf_error.SINGULAR, NULL)
        return zpack(nan, nan)
    elif z.real > SMALLX or fabs(z.imag) > SMALLY:
        return loggamma_stirling(z)
    elif zabs(z - 1) <= 0.5:
        return loggamma_taylor(z)
    elif zabs(z - 2) < 0.5:
        # Recurrence relation and the Taylor series around 1
        return zlog1(z - 1) + loggamma_taylor(z - 1)
    elif z.real < 0.1:
        # Reflection formula; see Proposition 3.1 in [1]
        tmp = npy_copysign(TWOPI, z.imag)*floor(0.5*z.real + 0.25)
        return zpack(LOGPI, tmp) - zlog(sinpi(z)) - loggamma(1 - z)
    elif npy_signbit(z.imag) == 0:
        # z.imag >= 0 but is not -0.0
        return loggamma_recurrence(z)
    else:
        return loggamma_recurrence(z.conjugate()).conjugate()
        

@cython.cdivision(True)
cdef inline double complex loggamma_recurrence(double complex z) nogil:
    """Backward recurrence relation; see Proposition 2.2 in [1]."""
    cdef:
        int n = <int>(SMALLX - z.real) + 1
        int signflips = 0
        int i
        double complex po, pn

    po = z + n - 1
    pn = po
    for i in range(2, n + 1):
        pn = po*(z + n - i)
        signflips += 1 if po.imag >=0 and pn.imag < 0 else 0
        po = pn
    return loggamma_stirling(z + n) - (zlog(pn) + signflips*TWOPI*1J)


@cython.cdivision(True)
cdef inline double complex loggamma_stirling(double complex z) nogil:
    cdef:
        int n = 1
        # Bernoulli numbers B_{2k} for 1 <= k <= 16
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
        double complex rzz = zdiv(coeff, z)
        double complex term

    res += coeff*bernoulli2k[n-1]/(2*n*(2*n - 1))
    for n in range(2, 17):
        coeff *= rzz
        term = coeff*bernoulli2k[n-1]/(2*n*(2*n - 1))
        res += term
        if zabs(term) <= TOL*zabs(res):
            break
    return res


@cython.cdivision(True)
cdef inline double complex loggamma_taylor(double complex z) nogil:
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
        if zabs(coeff/res) < TOL:
            break
    return res


cdef inline double complex cgamma(double complex z) nogil:
    """Compute Gamma(z) using loggamma."""
    if z.real <= 0 and z == floor(z.real):
        # Poles
        sf_error.error("gamma", sf_error.SINGULAR, NULL)
        return nan + 1J*nan
    return zexp(loggamma(z))


cdef inline double complex crgamma(double complex z) nogil:
    """Compute 1/Gamma(z) using loggamma."""
    if z.real <= 0 and z == floor(z.real):
        # Zeros at 0, -1, -2, ...
        return 0
    return zexp(-loggamma(z))
