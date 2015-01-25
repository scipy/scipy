
from libc.math cimport log, log1p, expm1, exp


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


cdef inline double inv_boxcox(double x, double lmbda) nogil:
    if lmbda == 0:
        return exp(x)
    else:
        return exp(log1p(lmbda * x) / lmbda)


cdef inline double inv_boxcox1p(double x, double lmbda) nogil:
    if lmbda == 0:
        return expm1(x)
    else:
        return expm1(log1p(lmbda * x) / lmbda)
