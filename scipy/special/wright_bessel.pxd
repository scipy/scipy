# -*-cython-*-
#
# Implementation of Wright's generalized Bessel function Phi [1].
# TODO: Based on the MPMath implementation [2], and documentation [3].
#
# Copyright: Christian Lorentzen
#
# Distributed under the same license as SciPy
#
# References:
# [1] https://dlmf.nist.gov/10.46.E1
# [2] mpmath source code, Subversion revision 990
#     http://code.google.com/p/mpmath/source/browse/trunk/mpmath/functions.py?spec=svn994&r=992
# [3] mpmath source code, Subversion revision 994
#     http://code.google.com/p/mpmath/source/browse/trunk/mpmath/function_docs.py?spec=svn994&r=994

# TODO:

import cython
from libc.math cimport exp

cdef extern from "_c99compat.h":
    int sc_isnan(double x) nogil
    int sc_isinf(double x) nogil

from ._cephes cimport rgamma
from ._complexstuff cimport inf, nan
from . cimport sf_error


@cython.cdivision(True)
cdef inline double wright_bessel_scalar(double rho, double beta, double x) nogil:
    cdef:
        int k
        double xk_k, res

    if sc_isnan(rho) or sc_isnan(beta) or sc_isnan(x):
        return nan
    elif rho < 0 or beta < 0 or x < 0:
        sf_error.error("wright_bessel", sf_error.DOMAIN, NULL)
        return nan
    elif sc_isinf(x):
        return inf
    elif sc_isinf(rho) or sc_isinf(beta):
        return 0
    elif x == 0:
        return rgamma(beta)
    elif rho == 0:
        return exp(x) * rgamma(beta)
    elif x <= 1:
        # rgamma(178) = 3e-323, rgamma(179) = 0.
        if beta >= 179:
            return 0
        else:
            # 18 term Taylor Series => error mostly smaller 1e-15
            # term k=0
            res = rgamma(beta)
            xk_k = 1.  # x^k/k!
            # term k=1,2,...
            for k in range(1, 18):
                xk_k *= x/k
                res += xk_k * rgamma(rho*k + beta)

            return res
    else:
        sf_error.error("wright_bessel", sf_error.DOMAIN, NULL)
