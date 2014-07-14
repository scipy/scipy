# -*- cython -*-
"""
Many Scipy special functions originally cast silently double input
arguments to integers.

Here, we define such unsafe wrappers manually.

"""

cimport sf_error

from sph_harm cimport sph_harmonic

cdef extern from "cephes.h":
    double bdtrc(int k, int n, double p) nogil
    double bdtr(int k, int n, double p) nogil 
    double bdtri(int k, int n, double y) nogil
    double expn(int n, double x) nogil
    double hyp2f0(double a, double b, double x, int type, double *err) nogil
    double nbdtrc(int k, int n, double p) nogil
    double nbdtr(int k, int n, double p) nogil 
    double nbdtri(int k, int n, double p) nogil
    double pdtrc(int k, double m) nogil
    double pdtr(int k, double m) nogil
    double pdtri(int k, double y) nogil
    double kn(int n, double x) nogil
    double yn(int n, double x) nogil
    double smirnov(int n, double e) nogil
    double smirnovi(int n, double p) nogil

cdef extern from "amos_wrappers.h":
    double cbesk_wrap_real_int(int n, double z) nogil

cdef extern from "Python.h":
    # Purposefully ignore the raised PyError --- assume the ufunc will collect it
    int PyErr_WarnEx_noerr "PyErr_WarnEx" (object, char *, int)

cdef inline void _legacy_cast_check(char *func_name, double x, double y) nogil:
    if <int>x != x or <int>y != y:
        with gil:
            PyErr_WarnEx_noerr(RuntimeWarning,
                               "floating point number truncated to an integer",
                               1)

cdef inline double complex sph_harmonic_unsafe(double m, double n, double theta, double phi) nogil:
    _legacy_cast_check("sph_harm", m, n)
    return sph_harmonic(<int>m, <int> n, theta, phi)
cdef inline double bdtrc_unsafe(double k, double n, double p) nogil:
    _legacy_cast_check("bdtrc", k, n)
    return bdtrc(<int>k, <int>n, p)
cdef inline double bdtr_unsafe(double k, double n, double p) nogil:
    _legacy_cast_check("bdtr", k, n)
    return bdtr(<int>k, <int>n, p)
cdef inline double bdtri_unsafe(double k, double n, double y) nogil:
    _legacy_cast_check("bdtri", k, n)
    return bdtri(<int>k, <int>n, y)
cdef inline double expn_unsafe(double n, double x) nogil:
    _legacy_cast_check("expn", n, 0)
    return expn(<int>n, x)
cdef inline double hyp2f0_unsafe(double a, double b, double x, double type, double *err) nogil:
    _legacy_cast_check("hyp2f0", type, 0)
    return hyp2f0(a, b, x, <int>type, err)
cdef inline double nbdtrc_unsafe(double k, double n, double p) nogil:
    _legacy_cast_check("nbdtrc", k, n)
    return nbdtrc(<int>k, <int>n, p)
cdef inline double nbdtr_unsafe(double k, double n, double p)  nogil:
    _legacy_cast_check("nbdtr", k, n)
    return nbdtr(<int>k, <int>n, p)
cdef inline double nbdtri_unsafe(double k, double n, double p) nogil:
    _legacy_cast_check("nbdtri", k, n)
    return nbdtri(<int>k, <int>n, p)
cdef inline double pdtrc_unsafe(double k, double m) nogil:
    _legacy_cast_check("pdtrc", k, 0)
    return pdtrc(<int>k, m)
cdef inline double pdtr_unsafe(double k, double m) nogil: 
    _legacy_cast_check("pdtr", k, 0)
    return pdtr(<int>k, m)
cdef inline double pdtri_unsafe(double k, double y) nogil:
    _legacy_cast_check("pdtri", k, 0)
    return pdtri(<int>k, y)
cdef inline double kn_unsafe(double n, double x) nogil:
    _legacy_cast_check("kn", n, 0)
    return cbesk_wrap_real_int(<int>n, x)
cdef inline double yn_unsafe(double n, double x) nogil:
    _legacy_cast_check("yn", n, 0)
    return yn(<int>n, x)
cdef inline double smirnov_unsafe(double n, double e) nogil:
    _legacy_cast_check("smirnov", n, 0)
    return smirnov(<int>n, e)
cdef inline double smirnovi_unsafe(double n, double p) nogil:
    _legacy_cast_check("smirnovi", n, 0)
    return smirnovi(<int>n, p)
