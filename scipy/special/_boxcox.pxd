
from libc.math cimport log, log1p, expm1


cdef inline double boxcox(double x, double lmbda) nogil:
    if lmbda == 0:
        return log(x)
    else:
        return expm1(lmbda * log(x)) / lmbda


cdef inline double boxcox1p(double x, double lmbda) nogil:
    if lmbda == 0:
        return log1p(x)
    else:
        return expm1(lmbda * log1p(x)) / lmbda
