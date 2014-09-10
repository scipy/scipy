cdef extern from "numpy/npy_math.h":
    double inf "NPY_INFINITY"

from libc.math cimport log

cdef inline double entr(double x) nogil:
    if x > 0:
        return -x * log(x)
    elif x == 0:
        return 0
    else:
        return -inf
