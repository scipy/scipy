# -*- cython -*-
"""
Many Scipy special functions originally cast silently double input
arguments to integers.

Here, we define such unsafe wrappers manually.

"""

from libc.math cimport isnan, isinf, NAN

from . cimport sf_error
from ._ellip_harm cimport ellip_harmonic
from .sph_harm cimport sph_harmonic
from ._cephes cimport (bdtrc, bdtr, bdtri, expn, nbdtrc,
                       nbdtr, nbdtri, pdtri, kn, yn,
                       smirnov, smirnovi, smirnovc, smirnovci, smirnovp)

cdef extern from "amos_wrappers.h":
    double cbesk_wrap_real_int(int n, double z) nogil

cdef extern from "Python.h":
    # Purposefully ignore the raised PyError --- assume the ufunc will collect it
    int PyErr_WarnEx_noerr "PyErr_WarnEx" (object, char *, int)

cdef inline void _legacy_cast_check(char *func_name, double x, double y) noexcept nogil:
    if <int>x != x or <int>y != y:
        with gil:
            PyErr_WarnEx_noerr(RuntimeWarning,
                               "floating point number truncated to an integer",
                               1)

cdef inline void _legacy_deprecation(char *func_name, double x, double y) noexcept nogil:
        with gil:
            PyErr_WarnEx_noerr(DeprecationWarning,
                               "non-integer arg n is deprecated, removed in SciPy 1.7.x",
                               1)

cdef inline double complex sph_harmonic_unsafe(double m, double n,
                                               double theta, double phi) noexcept nogil:
    if isnan(m) or isnan(n):
        return NAN
    _legacy_cast_check("sph_harm", m, n)
    return sph_harmonic(<int>m, <int> n, theta, phi)

cdef inline double ellip_harmonic_unsafe(double h2, double k2, double n,
                                         double p, double l, double signm,
                                         double signn) noexcept nogil:
    if isnan(n) or isnan(p):
        return NAN
    _legacy_cast_check("_ellip_harm", n, p)
    return ellip_harmonic(h2, k2, <int>n, <int>p, l, signm, signn)

cdef inline double bdtr_unsafe(double k, double n, double p) noexcept nogil:
    _legacy_deprecation("bdtr", k, n)
    if isnan(n) or isinf(n):
        return NAN
    else:
        return bdtr(k, <int>n, p)

cdef inline double bdtrc_unsafe(double k, double n, double p) noexcept nogil:
    _legacy_deprecation("bdtrc", k, n)
    if isnan(n) or isinf(n):
        return NAN
    else:
        return bdtrc(k, <int>n, p)

cdef inline double bdtri_unsafe(double k, double n, double p) noexcept nogil:
    _legacy_deprecation("bdtri", k, n)
    if isnan(n) or isinf(n):
        return NAN
    else:
        return bdtri(k, <int>n, p)

cdef inline double expn_unsafe(double n, double x) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("expn", n, 0)
    return expn(<int>n, x)

cdef inline double nbdtrc_unsafe(double k, double n, double p) noexcept nogil:
    if isnan(k) or isnan(n):
        return NAN
    _legacy_cast_check("nbdtrc", k, n)
    return nbdtrc(<int>k, <int>n, p)

cdef inline double nbdtr_unsafe(double k, double n, double p) noexcept nogil:
    if isnan(k) or isnan(n):
        return NAN
    _legacy_cast_check("nbdtr", k, n)
    return nbdtr(<int>k, <int>n, p)

cdef inline double nbdtri_unsafe(double k, double n, double p) noexcept nogil:
    if isnan(k) or isnan(n):
        return NAN
    _legacy_cast_check("nbdtri", k, n)
    return nbdtri(<int>k, <int>n, p)

cdef inline double pdtri_unsafe(double k, double y) noexcept nogil:
    if isnan(k):
        return k
    _legacy_cast_check("pdtri", k, 0)
    return pdtri(<int>k, y)

cdef inline double kn_unsafe(double n, double x) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("kn", n, 0)
    return cbesk_wrap_real_int(<int>n, x)

cdef inline double yn_unsafe(double n, double x) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("yn", n, 0)
    return yn(<int>n, x)

cdef inline double smirnov_unsafe(double n, double e) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("smirnov", n, 0)
    return smirnov(<int>n, e)

cdef inline double smirnovc_unsafe(double n, double e) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("smirnovc", n, 0)
    return smirnovc(<int>n, e)

cdef inline double smirnovp_unsafe(double n, double e) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("smirnovp", n, 0)
    return smirnovp(<int>n, e)

cdef inline double smirnovi_unsafe(double n, double p) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("smirnovi", n, 0)
    return smirnovi(<int>n, p)

cdef inline double smirnovci_unsafe(double n, double p) noexcept nogil:
    if isnan(n):
        return n
    _legacy_cast_check("smirnovci", n, 0)
    return smirnovci(<int>n, p)
