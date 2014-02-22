
cdef extern from "numpy/npy_math.h":
    double inf "NPY_INFINITY"

from libc.math cimport log


cdef inline double kl_div(double x, double y) nogil:
    if x < 0 or y < 0 or (y == 0 and x != 0):
        # extension of the natural domain to preserve convexity
        return inf;
    elif x == inf or y == inf:
        # limits within the natural domain
        return inf;
    elif x == 0:
        return y;
    else:
        return x*log(x) - x*log(y) - x + y;
