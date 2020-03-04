import ctypes
from libc.math cimport sqrt, fabs
from libc.stdlib cimport free
from numpy import nan

cdef extern from "Python.h":
    object PyCapsule_New(void *pointer, char *name, void *destructor)

from scipy._lib._ccallback import LowLevelCallable
from ._ellip_harm cimport ellip_harmonic, ellip_harm_eval, lame_coefficients

ctypedef struct _ellip_data_t:
    double *eval
    double h2, k2
    int n, p

cdef double _F_integrand(double t, void *user_data) nogil:
    cdef _ellip_data_t *data = <_ellip_data_t *>user_data
    cdef double h2, k2, t2, i, a, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = data[0].h2
    k2 = data[0].k2
    n = data[0].n
    p = data[0].p
    eval = data[0].eval
    i = ellip_harm_eval( h2, k2, n, p, 1/t, eval, 1, 1)
    result = 1/(i*i*sqrt(1 - t2*k2)*sqrt(1 - t2*h2))
    return result

cdef double _F_integrand1(double t, void *user_data) nogil:
    cdef _ellip_data_t *data = <_ellip_data_t *>user_data
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = data[0].h2
    k2 =data[0].k2
    n = data[0].n
    p = data[0].p
    eval = data[0].eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = i*i/sqrt((t + h)*(t + k))
    return result

cdef double _F_integrand2(double t, void *user_data) nogil:
    cdef _ellip_data_t *data = <_ellip_data_t *>user_data
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = data[0].h2
    k2 =data[0].k2
    n = data[0].n
    p = data[0].p
    eval = data[0].eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = t2*i*i/sqrt((t + h)*(t + k))
    return result

cdef double _F_integrand3(double t, void *user_data) nogil:
    cdef _ellip_data_t *data = <_ellip_data_t *>user_data
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double * eval
    t2 = t*t
    h2 = data[0].h2
    k2 =data[0].k2
    n = data[0].n
    p = data[0].p
    eval = data[0].eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = i*i/sqrt((t + h)*(k2 - t2))
    return result

cdef double _F_integrand4(double t, void *user_data) nogil:
    cdef _ellip_data_t *data = <_ellip_data_t *>user_data
    cdef double h2, k2, t2, i, a, h, result
    cdef int n, p
    cdef double *eval
    t2 = t*t
    h2 = data[0].h2
    k2 =data[0].k2
    n = data[0].n
    p = data[0].p
    eval = data[0].eval

    h = sqrt(h2)
    k = sqrt(k2)
    i = ellip_harm_eval( h2, k2, n, p, t, eval, 1, 1)
    result = i*i*t2/sqrt((t + h)*(k2 - t2))
    return result


def _ellipsoid(double h2, double k2, int n, int p, double s):
    import scipy.special._ellip_harm_2 as mod
    from scipy.integrate import quad

    cdef _ellip_data_t data

    cdef double * eval
    cdef void *bufferp
    eval = lame_coefficients(h2, k2, n, p, &bufferp, 1, 1)
    if not eval:
        return nan

    data.h2 = h2
    data.k2 = k2
    data.n = n
    data.p = p
    data.eval = eval

    cdef double res, err

    try:
        capsule = PyCapsule_New(<void*>&data, NULL, NULL)
        res, err = quad(LowLevelCallable.from_cython(mod, "_F_integrand", capsule), 0, 1/s,
                                                     epsabs=1e-300, epsrel=1e-15)
    finally:
        free(bufferp)
    if err > 1e-10*fabs(res) + 1e-290:
        return nan
    res = res*(2*n + 1)*ellip_harmonic( h2, k2, n, p, s, 1, 1)
    return res


def _ellipsoid_norm(double h2, double k2, int n, int p):
    import scipy.special._ellip_harm_2 as mod
    from scipy.integrate import quad

    cdef _ellip_data_t data

    cdef double *eigv
    cdef void *bufferp
    eval = lame_coefficients(h2, k2, n, p, &bufferp, 1, 1)
    if not eval:
        return nan

    data.h2 = h2
    data.k2 = k2
    data.n = n
    data.p = p
    data.eval = eval

    cdef double res, res1, res2, res3, err, err1, err2, err3

    h = sqrt(h2)
    k = sqrt(k2)
    try:
        capsule = PyCapsule_New(<void*>&data, NULL, NULL)

        wvar = (-0.5, -0.5)

        res, err = quad(LowLevelCallable.from_cython(mod, "_F_integrand1", capsule), h, k,
                        epsabs=1e-300, epsrel=1e-15, weight="alg", wvar=wvar)

        res1, err1 = quad(LowLevelCallable.from_cython(mod, "_F_integrand2", capsule), h, k,
                          epsabs=1e-300, epsrel=1e-15, weight="alg", wvar=wvar)

        wvar = (0, -0.5)

        res2, err2 = quad(LowLevelCallable.from_cython(mod, "_F_integrand3", capsule), 0, h,
                          epsabs=1e-300, epsrel=1e-15, weight="alg", wvar=wvar)

        res3, err3 = quad(LowLevelCallable.from_cython(mod, "_F_integrand4", capsule), 0, h,
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
