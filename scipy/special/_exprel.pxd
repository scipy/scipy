#cython: cdivison=True
cdef extern from "numpy/npy_math.h":
    double NPY_INFINITY

from libc.math cimport expm1, fabs, log

cdef inline double exprel(double x) nogil:
    if fabs(x) < 1e-16:
        return 1.0
    elif x > 717:  # near log(DBL_MAX)
        return NPY_INFINITY
    else:
        return expm1(x) / x

