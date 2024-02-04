from libc.math cimport log1p

from ._complexstuff cimport zlog, zisnan, number_t
from ._cunity cimport clog1p

cdef inline number_t xlogy(number_t x, number_t y) noexcept nogil:
    if x == 0 and not zisnan(y):
        return 0
    else:
        return x * zlog(y)

cdef inline number_t xlog1py(number_t x, number_t y) noexcept nogil:
    if x == 0 and not zisnan(y):
        return 0
    else:
        if number_t is double:
            return x * log1p(y)
        else:
            return x * clog1p(y)
