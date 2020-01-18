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
from libc.math cimport exp, floor

cdef extern from "_c99compat.h":
    int sc_isnan(double x) nogil
    int sc_isinf(double x) nogil

from ._cephes cimport Gamma, rgamma
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
    elif x <= 1:
        # series expansion until term k such that a*k+b <= rgamma_zero
        k_max = int(floor((rgamma_zero - b) / a))

        # 19 term Taylor Series => error mostly smaller 1e-15
        if k_max > 19:
            k_max = 19
        return _wb_series(a, b, x, 0, k_max)
    else:
        sf_error.error("wright_bessel", sf_error.DOMAIN, NULL)
