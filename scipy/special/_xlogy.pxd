# -*-cython-*-

from _complexstuff cimport zlog, npy_log1p, zisnan, number_t
from _cunity cimport clog1p

cdef inline number_t xlogy(number_t x, number_t y) nogil:
    if x == 0 and not zisnan(y):
        return 0
    else:
        return x * zlog(y)

cdef inline number_t xlog1py(number_t x, number_t y) nogil:
    if x == 0 and not zisnan(y):
        return 0
    else:
        if number_t is double:
            return x * npy_log1p(y)
        else:
            return x * clog1p(y)
