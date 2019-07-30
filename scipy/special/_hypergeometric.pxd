from . cimport sf_error

from ._cephes cimport poch

cdef extern from "numpy/npy_math.h":
    double NPY_NAN
    double NPY_INFINITY

cdef extern from 'specfun_wrappers.h':
    double hypU_wrap(double, double, double) nogil


cdef inline double hyperu(double a, double b, double x) nogil:
    if x < 0.0:
        sf_error.error("hyperu", sf_error.DOMAIN, NULL)
        return NPY_NAN

    if x == 0.0:
        if b > 1.0:
            # DMLF 13.2.16-18
            sf_error.error("hyperu", sf_error.SINGULAR, NULL)
            return NPY_INFINITY
        else:
            # DLMF 13.2.14-15 and 13.2.19-21
            return poch(1.0 - b + a, -a)

    return hypU_wrap(a, b, x)
