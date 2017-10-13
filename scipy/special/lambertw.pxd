# -*-cython-*-
#
# Implementation of the Lambert W function [1]. Based on the MPMath
# implementation [2], and documentaion [3].
#
# Copyright: Yosef Meller, 2009
# Author email: mellerf@netvision.net.il
#
# Distributed under the same license as SciPy
#
# References:
# [1] On the Lambert W function, Adv. Comp. Math. 5 (1996) 329-359,
#     available online: http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf
# [2] mpmath source code, Subversion revision 990
#     http://code.google.com/p/mpmath/source/browse/trunk/mpmath/functions.py?spec=svn994&r=992
# [3] mpmath source code, Subversion revision 994
#     http://code.google.com/p/mpmath/source/browse/trunk/mpmath/function_docs.py?spec=svn994&r=994

# TODO: use a series expansion when extremely close to the branch point
# at `-1/e` and make sure that the proper branch is chosen there

import cython

from . cimport sf_error

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil

from ._complexstuff cimport *

DEF EXPN1 = 0.36787944117144232159553  # exp(-1)


@cython.cdivision(True)
cdef inline double complex lambertw_scalar(double complex z, long k, double tol) nogil:
    cdef int i
    cdef double absz
    cdef double complex w
    cdef double complex ew, wew, wewz, wn

    if zisnan(z):
        return z
    elif z.real == inf:
        return z + 2*k*pi*1j
    elif z.real == -inf:
        return -z + (2*k+1)*pi*1j
    elif z == 0:
        if k == 0:
            return z
        sf_error.error("lambertw", sf_error.SINGULAR, NULL)
        return -inf

    absz = zabs(z)
    # Get an initial guess for Halley's method
    if k == 0:
        if absz <= EXPN1:
            w = z
        elif z.imag and absz <= 0.7:
            if zabs(z + 0.5) < 0.1:
                if z.imag > 0:
                    w = 0.7 + 0.7j
                else:
                    w = 0.7 - 0.7j
            else:
                w = z
        else:
            w = zlog(z) + 2*pi*k*1j
    elif k == -1:
        if absz <= EXPN1 and z.imag == 0 and z.real < 0:
            w = log(-z.real)
        else:
            w = zlog(z) + 2*pi*k*1j
    else:
        w = zlog(z) + 2*pi*k*1j

    # Halley's method
    for i in range(100):
        ew = zexp(w)
        wew = w*ew
        wewz = wew - z
        wn = w - wewz / (wew + ew - (w + 2)*wewz/(2*w + 2))
        if zabs(wn - w) < tol*zabs(wn):
            return wn
        else:
            w = wn

    sf_error.error("lambertw", sf_error.SLOW,
                   "iteration failed to converge: %g + %gj",
                   <double>z.real, <double>z.imag)
    return nan
