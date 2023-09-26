# -*-cython-*-
#
# Implementation of the Lambert W function [1]. Based on the MPMath
# implementation [2], and documentation [3].
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
from libc.math cimport exp, log, INFINITY, NAN, M_PI, M_E

from . cimport sf_error
from ._evalpoly cimport cevalpoly
from ._complexstuff cimport *

DEF EXPN1 = 0.36787944117144232159553  # exp(-1)


@cython.cdivision(True)
cdef inline double complex lambertw_scalar(double complex z, long k, double tol) noexcept nogil:
    cdef int i
    cdef double absz, p
    cdef double complex w
    cdef double complex ew, wew, wewz, wn
    cdef double OMEGA = 0.56714329040978387299997  # W(1, 0)

    if zisnan(z):
        return z
    elif z.real == INFINITY:
        return z + 2*M_PI*k*1j
    elif z.real == -INFINITY:
        return -z + (2*M_PI*k+M_PI)*1j
    elif z == 0:
        if k == 0:
            return z
        sf_error.error("lambertw", sf_error.SINGULAR, NULL)
        return -INFINITY
    elif z == 1 and k == 0:
        # Split out this case because the asymptotic series blows up
        return OMEGA

    absz = zabs(z)
    # Get an initial guess for Halley's method
    if k == 0:
        if zabs(z + EXPN1) < 0.3:
            w = lambertw_branchpt(z)
        elif (-1.0 < z.real < 1.5 and zabs(z.imag) < 1.0
              and -2.5*zabs(z.imag) - 0.2 < z.real):
            # Empirically determined decision boundary where the Pade
            # approximation is more accurate.
            w = lambertw_pade0(z)
        else:
            w = lambertw_asy(z, k)
    elif k == -1:
        if absz <= EXPN1 and z.imag == 0 and z.real < 0:
            w = log(-z.real)
        else:
            w = lambertw_asy(z, k)
    else:
        w = lambertw_asy(z, k)

    # Halley's method; see 5.9 in [1]
    if w.real >= 0:
        # Rearrange the formula to avoid overflow in exp
        for i in range(100):
            ew = zexp(-w)
            wewz = w - z*ew
            wn = w - wewz/(w + 1 - (w + 2)*wewz/(2*w + 2))
            if zabs(wn - w) <= tol*zabs(wn):
                return wn
            else:
                w = wn
    else:
        for i in range(100):
            ew = zexp(w)
            wew = w*ew
            wewz = wew - z
            wn = w - wewz/(wew + ew - (w + 2)*wewz/(2*w + 2))
            if zabs(wn - w) <= tol*zabs(wn):
                return wn
            else:
                w = wn

    sf_error.error("lambertw", sf_error.SLOW,
                   "iteration failed to converge: %g + %gj",
                   <double>z.real, <double>z.imag)
    return zpack(NAN, NAN)


@cython.cdivision(True)
cdef inline double complex lambertw_branchpt(double complex z) noexcept nogil:
    """Series for W(z, 0) around the branch point; see 4.22 in [1]."""
    cdef double *coeffs = [-1.0/3.0, 1.0, -1.0]
    cdef double complex p = zsqrt(2*(M_E*z + 1))

    return cevalpoly(coeffs, 2, p)


@cython.cdivision(True)
cdef inline double complex lambertw_pade0(double complex z) noexcept nogil:
    """(3, 2) Pade approximation for W(z, 0) around 0."""
    cdef:
        double *num = [
            12.85106382978723404255,
            12.34042553191489361902,
            1.0
        ]
        double *denom = [
            32.53191489361702127660,
            14.34042553191489361702,
            1.0
        ]

    # This only gets evaluated close to 0, so we don't need a more
    # careful algorithm that avoids overflow in the numerator for
    # large z.
    return z*cevalpoly(num, 2, z)/cevalpoly(denom, 2, z)


@cython.cdivision(True)
cdef inline double complex lambertw_asy(double complex z, long k) noexcept nogil:
    """Compute the W function using the first two terms of the
    asymptotic series. See 4.20 in [1].

    """
    cdef double complex w

    w = zlog(z) + 2*M_PI*k*1j
    return w - zlog(w)
