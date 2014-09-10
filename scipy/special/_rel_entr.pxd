cdef extern from "numpy/npy_math.h":
    double inf "NPY_INFINITY"

from libc.math cimport log

cdef inline double rel_entr(double x, double y) nogil:
    if x > 0 and y > 0:
        return x * log(x / y)
    elif x == 0 and y >= 0:
        return 0
    else:
        return inf
