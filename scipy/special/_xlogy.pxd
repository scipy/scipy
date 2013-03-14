from _complexstuff cimport zlog, zisnan, number_t

cdef inline number_t xlogy(number_t y, number_t x) nogil:
    if y == 0 and not zisnan(x):
        return 0
    else:
        return y * zlog(x)
