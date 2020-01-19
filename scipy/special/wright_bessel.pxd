# -*-cython-*-
#
# Implementation of Wright's generalized Bessel function Phi [1].
# TODO: Based on the MPMath implementation [2].
#
# Copyright: Christian Lorentzen
#
# Distributed under the same license as SciPy
#
# Implementation Overview:
# - Implement only for non-negative values of rho=a, beta=b and z=x
# - x <= 1: Use series expansion [1]
# - large a or b: Use series expansion
# - large x: Use asymptotic expansion [3, 4]
#     Gets worse if rho=a gets large
# - intermediate x (all the rest): integral representation
#
# References:
# [1] https://dlmf.nist.gov/10.46.E1
# [2] TODO: mpmath source code or notebook
# [3] E. M. Wright (1935), The asymptotic expansion of the generalized Bessel
#     function. Proc. London Math. Soc. (2) 38, pp. 257–270.
#     https://doi.org/10.1112/plms/s2-38.1.257
# [4] R. B. Paris (2017), The asymptotics of the generalised Bessel function,
#     Mathematica Aeterna, Vol. 7, 2017, no. 4, 381 - 406,
#     https://arxiv.org/abs/1711.03006
# [5] Y. F. Luchko (2008), Algorithms for Evaluation of the Wright Function for
#     the Real Arguments’ Values, Fractional Calculus and Applied Analysis 11(1)
#     http://sci-gems.math.bas.bg/jspui/bitstream/10525/1298/1/fcaa-vol11-num1-2008-57p-75p.pdf

import cython
from libc.math cimport exp, floor, pow, sqrt, M_PI

cdef extern from "_c99compat.h":
    int sc_isnan(double x) nogil
    int sc_isinf(double x) nogil

from ._cephes cimport Gamma, rgamma, zeta
from ._digamma cimport digamma
from ._complexstuff cimport inf, nan
from . cimport sf_error

# rgamma_zero: smallest value x for which rgamma(x) == 0 as x gets large
DEF rgamma_zero = 178.47241115886637


@cython.cdivision(True)
cdef inline double _wb_series(double a, double b, double x,
    unsigned int nstart, unsigned int nstop) nogil:
    cdef:
        unsigned int k
        double xk_k, res

    xk_k = x**nstart * Gamma(nstart+1)  # x^k/k!
    res = xk_k * rgamma(nstart*a + b)
    # term k=nstart+1 , +2 +3, ...
    if nstop > nstart:
        for k in range(nstart+1, nstop):
            xk_k *= x/k
            res += xk_k * rgamma(a*k + b)

    return res


@cython.cdivision(True)
cdef inline double _wb_small_a(double a, double b, double x) nogil:
    """Compute series of Wright bessel function in 'a' to order 2."""
    cdef double dg, dg1, res
    dg = digamma(b)
    # dg1 = polygamma(1, b)
    dg1 = zeta(2, x)
    # res = 1 - a*x*dg + a**2/2*x*(1+x)*(dg**2 - dg1)
    res = 1 + a*x*(-dg + 0.5*a*(1+x)*(dg**2 - dg1))
    res *= exp(x) * rgamma(b)
    return res


@cython.cdivision(True)
cdef inline double _wb_asymptotic(double a, double b, double x) nogil:
    """Compute asymptotic series for Wright Bessel function.

    Asymptotic series for positive a according to Wright (1935).
    Z = (a*x)**(1/(1+a))
    Phi ~ Z**(1/2-b) * exp((1+a)/a * Z) * sum_k (-1)**k * a_k / Z**k
    For coefficients up to a_3, see Paris (2017).
    """
    cdef:
        double A[7]  # powers of a
        double B[9]  # powers of b
        double Ap1[5]  # powers of (1+a)
        double C[5]  # coefficients of asymptotic series a_k
        double Z, Zp, res
        int k

    A[0] = 1.
    B[0] = 1.
    Ap1[0] = 1.
    for k in range(1, 7):
        A[k] = A[k-1] * a
    for k in range(1, 9):
        B[k] = B[k-1] * b
    for k in range(1, 5):
        Ap1[k] = Ap1[k-1] * (1 + a)

    C[0] = 1./sqrt(2. * M_PI * Ap1[1])

    C[1] = C[0] / (24 * Ap1[1])
    C[1] *= (2*a + 1)*(2 + a) - 12*b*(1 + a - b)

    C[2] = C[0] / (1152 * Ap1[2])
    C[2] *= ((2 + a)*(1 + 2*a)*(2 - 19*a + 2*A[2]) \
            + 24*b*Ap1[1]*(2 + 7*a - 6*A[2]) \
            - 24*B[2]*(4 - 5*a - 20*A[2]) - 96*B[3]*(1 + 5*a) + 144*B[4])

    C[3] = -C[0] / (414720*Ap1[3])
    C[3] *= (2 + a)*(1 + 2*a)*(556 + 1628*a - 9093*A[2] + 1628*A[3] \
                               + 556*A[4]) \
            - 180*b*Ap1[1]*(12 - 172*a - 417*A[2] + 516*A[3] - 20*A[4]) \
            - 180*B[2]*(76 + 392*a - 567*A[2] - 1288*A[3] + 364*A[4]) \
            + 1440*B[3]*(8 - 63*a - 147*A[2] + 112*A[3]) \
            + 10800*B[4]*(2 + 7*a - 14*A[2]) \
            - 8640*B[5]*(1 - 7*a) - 8640*B[6]

    C[4] = C[0] / (39813120*Ap1[4])
    C[4] *= 103680*B[8] - 414720*B[7]*(3*a - 1) + 725760*B[6]*a*(8*a - 7) \
            - 48384*B[5]*(274*A[3] - 489*A[2] + 39*a + 26) \
            + 30240*B[4]*(500*A[4] - 1740*A[3] + 495*A[2] + 340*a - 12) \
            - 2880*B[3]*(2588*A[5] - 19780*A[4] + 14453*A[3] + 9697*A[2] \
                         - 1892*a - 404) \
            + 48*B[2]*(11488*A[6] - 547836*A[5] + 1007484*A[4] + 593353*A[3] \
                       - 411276*A[2] - 114396*a + 4288) \
            + 48*b*Ap1[1]*(7784*A[6] + 48180*A[5] - 491202*A[4] + 336347*A[3] \
                           + 163734*A[2] - 28908*a - 5560) \
            - (a + 2)*(2*a + 1)*(4568*A[6] - 226668*A[5] - 465702*A[4] \
                                 + 2013479*A[3] - 465702*A[2] - 226668*a \
                                 + 4568)

    Z = pow(a * x, 1/Ap1[1])
    Zp = 1.
    res = C[0]
    for k in range(1, 5):
        Zp /= Z
        res += (-1)**k * C[k] * Zp
    res *= pow(Z, 0.5 - b) * exp(Ap1[1]/a * Z)


@cython.cdivision(True)
cdef inline double wright_bessel_scalar(double a, double b, double x) nogil:
    cdef:
        double xk_k, res
        unsigned int k_max

    if sc_isnan(a) or sc_isnan(b) or sc_isnan(x):
        return nan
    elif a < 0 or b < 0 or x < 0:
        sf_error.error("wright_bessel", sf_error.DOMAIN, NULL)
        return nan
    elif sc_isinf(x):
        return inf
    elif sc_isinf(a) or sc_isinf(b):
        return 0
    elif x == 0:
        return rgamma(b)
    elif a == 0:
        return exp(x) * rgamma(b)
    elif a <= 1e-5 and x <= 10:
        # series expansion in 'a' to order 2 has precision ~ 2e-12
        return _wb_small_a(a, b, x)
    elif x <= 1:
        # series expansion until term k such that a*k+b <= rgamma_zero
        k_max = int(floor((rgamma_zero - b) / a))

        # 18 term Taylor Series => error mostly smaller 1e-15
        if k_max > 18:
            k_max = 18
        return _wb_series(a, b, x, 0, k_max)
    elif x <= 2:
        # series expansion until term k such that a*k+b <= rgamma_zero
        k_max = int(floor((rgamma_zero - b) / a))

        # 20 term Taylor Series => error mostly smaller 1e-12 to 1e-13
        if k_max > 20:
            k_max = 20
        return _wb_series(a, b, x, 0, k_max)
    else:
        sf_error.error("wright_bessel", sf_error.DOMAIN, NULL)
