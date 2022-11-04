from libc.math cimport expm1, fabs, log, INFINITY

cdef inline double exprel(double x) nogil:
    if fabs(x) < 1e-16:
        return 1.0
    elif x > 717:  # near log(DBL_MAX)
        return INFINITY
    else:
        return expm1(x) / x

