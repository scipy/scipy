import ctypes
from libc.math cimport sqrt, fabs
from libc.stdlib cimport free
from numpy import nan
import scipy.integrate

# The access to global variables is protected by
# is protected by _ellip_lock in _ellip_harm.py

cdef double * _global_eval
cdef double _global_h2, _global_k2
cdef int _global_n, _global_p

from ._ellip_harm cimport ellip_harmonic, ellip_harm_eval, lame_coefficients

cdef double _F_integrand(double t) nogil:
    cdef double h2, k2, t2, i, a, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = _global_h2
    k2 = _global_k2
    n = _global_n
    p = _global_p
    eval = _global_eval
    i = ellip_harm_eval( h2, k2, n, p, 1/t, eval, 1, 1)
    result = 1/(i*i*sqrt(1 - t2*k2)*sqrt(1 - t2*h2))
    return result

_F_integrand_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
_F_integrand_ctypes = ctypes.cast(<size_t>&_F_integrand, _F_integrand_t)

cdef double _F_integrand1(double t) nogil:
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = _global_h2
    k2 =_global_k2
    n = _global_n
    p = _global_p
    eval = _global_eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = i*i/sqrt((t + h)*(t + k))
    return result

_F_integrand1_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
_F_integrand1_ctypes = ctypes.cast(<size_t>&_F_integrand1, _F_integrand1_t)

cdef double _F_integrand2(double t) nogil:
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = _global_h2
    k2 =_global_k2
    n = _global_n
    p = _global_p
    eval = _global_eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = t2*i*i/sqrt((t + h)*(t + k))
    return result

_F_integrand2_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
_F_integrand2_ctypes = ctypes.cast(<size_t>&_F_integrand2, _F_integrand2_t)

cdef double _F_integrand3(double t) nogil:
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = _global_h2
    k2 =_global_k2
    n = _global_n
    p = _global_p
    eval = _global_eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = i*i/sqrt((t + h)*(k2 - t2))
    return result

_F_integrand3_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
_F_integrand3_ctypes = ctypes.cast(<size_t>&_F_integrand3, _F_integrand3_t)

cdef double _F_integrand4(double t) nogil:
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double *eval
    t2 = t*t
    h2 = _global_h2
    k2 =_global_k2
    n = _global_n
    p = _global_p
    eval = _global_eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = i*i*t2/sqrt((t + h)*(k2 - t2))
    return result

_F_integrand4_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
_F_integrand4_ctypes = ctypes.cast(<size_t>&_F_integrand4, _F_integrand4_t)


def _ellipsoid(double h2, double k2, int n, int p, double s):
    global _global_h2
    global _global_k2
    global _global_n
    global _global_p
    global _global_eval

    cdef double * eval
    cdef void *bufferp
    eval = lame_coefficients(h2, k2, n, p, &bufferp, 1, 1)
    if not eval:
        return nan

    _global_h2 = h2
    _global_k2 = k2
    _global_n = n
    _global_p = p
    _global_eval = eval

    cdef double res, err

    try:
        res, err = scipy.integrate.quad(_F_integrand_ctypes, 0, 1/s,
                                        epsabs=1e-300, epsrel=1e-15)
    finally:
        free(bufferp)
    if err > 1e-10*fabs(res) + 1e-290:
        return nan
    res = res*(2*n + 1)*ellip_harmonic( h2, k2, n, p, s, 1, 1)
    return res


def _ellipsoid_norm(double h2, double k2, int n, int p):
    global _global_h2
    global _global_k2
    global _global_n
    global _global_p
    global _global_eval

    cdef double *eigv
    cdef void *bufferp
    eval = lame_coefficients(h2, k2, n, p, &bufferp, 1, 1)
    if not eval:
        return nan

    _global_h2 = h2
    _global_k2 = k2
    _global_n = n
    _global_p = p
    _global_eval = eval

    cdef double res, res1, res2, res3, err, err1, err2, err3

    h = sqrt(h2)
    k = sqrt(k2)
    try:
        quad = scipy.integrate.quad

        wvar = (-0.5, -0.5)

        res, err = quad(_F_integrand1_ctypes, h, k,
                        epsabs=1e-300, epsrel=1e-15, weight="alg", wvar=wvar)

        res1, err1 = quad(_F_integrand2_ctypes, h, k,
                          epsabs=1e-300, epsrel=1e-15, weight="alg", wvar=wvar)

        wvar = (0, -0.5)

        res2, err2 = quad(_F_integrand3_ctypes, 0, h,
                          epsabs=1e-300, epsrel=1e-15, weight="alg", wvar=wvar)

        res3, err3 = quad(_F_integrand4_ctypes, 0, h,
                          epsabs=1e-300, epsrel=1e-15, weight="alg", wvar=wvar)

    finally:
        free(bufferp)

    error = 8*(res2*err1 + err2*res1 + res*err3 + res3*err)
    result = 8*(res1*res2 - res*res3)

    if error > 10e-8*fabs(result):
        return nan
    return result


# Needed for the _sf_error calls in _ellip_harm.pxd

cimport numpy as np

np.import_array()
np.import_ufunc()

cdef extern from "numpy/ufuncobject.h":
    int PyUFunc_getfperr() nogil

cdef public int wrap_PyUFunc_getfperr() nogil:
    """
    Call PyUFunc_getfperr in a context where PyUFunc_API array is initialized;
    this avoids messing with the UNIQUE_SYMBOL #defines
    """
    return PyUFunc_getfperr()
