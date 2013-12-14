
from libc.math cimport log, pow


cdef inline double boxcox(double x, double lmbda) nogil:
    if lmbda == 0:
        return log(x)
    else:
        return (pow(x, lmbda) - 1.0) / lmbda
