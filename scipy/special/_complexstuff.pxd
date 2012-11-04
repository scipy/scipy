# -*-cython-*-
#
# Common functions required when doing complex arithmetic with Cython.
#

cimport numpy as np

cdef extern from "_complexstuff.h":
    double npy_cabs(np.npy_cdouble z) nogil
    np.npy_cdouble npy_clog(np.npy_cdouble z) nogil
    np.npy_cdouble npy_cexp(np.npy_cdouble z) nogil
    int npy_isnan(double x) nogil
    double inf "NPY_INFINITY"
    double pi "NPY_PI"
    double nan "NPY_NAN"

cdef inline bint zisnan(double complex x) nogil:
    return npy_isnan(x.real) or npy_isnan(x.imag)

cdef inline double zabs(double complex x) nogil:
    cdef double r
    r = npy_cabs((<np.npy_cdouble*>&x)[0])
    return r

cdef inline double complex zlog(double complex x) nogil:
    cdef np.npy_cdouble r
    r = npy_clog((<np.npy_cdouble*>&x)[0])
    return (<double complex*>&r)[0]

cdef inline double complex zexp(double complex x) nogil:
    cdef np.npy_cdouble r
    r = npy_cexp((<np.npy_cdouble*>&x)[0])
    return (<double complex*>&r)[0]
