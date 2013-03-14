# -*-cython-*-

from _complexstuff cimport zlog, npy_log1p, zisnan, number_t

cdef inline number_t xlogy(number_t x, number_t y) nogil:
    if x == 0 and not zisnan(y):
        return 0
    else:
        return x * zlog(y)

cdef inline double xlog1py(double x, double y) nogil:
    if x == 0 and not zisnan(y):
        return 0
    else:
        return x * npy_log1p(y)
