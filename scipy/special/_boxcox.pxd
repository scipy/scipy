
from libc.math cimport log, log1p, expm1, exp, fabs


cdef inline double boxcox(double x, double lmbda) nogil:
    # if lmbda << 1 and log(x) < 1.0, the lmbda*log(x) product can lose
    # precision, furthermore, expm1(x) == x for x < eps.
    # For doubles, the range of log is -744.44 to +709.78, with eps being
    # the smallest value produced.  This range means that we will have
    # abs(lmbda)*log(x) < eps whenever abs(lmbda) <= eps/-log(min double)
    # which is ~2.98e-19.  
    if fabs(lmbda) < 1e-19:
        return log(x)
    else:
        return expm1(lmbda * log(x)) / lmbda


cdef inline double boxcox1p(double x, double lmbda) nogil:
    # The argument given above in boxcox applies here with the modification
    # that the smallest value produced by log1p is the minimum representable
    # value, rather than eps.  The second condition here prevents unflow
    # when log1p(x) is < eps.
    cdef double lgx = log1p(x)
    if fabs(lmbda) < 1e-19 or (fabs(lgx) < 1e-289 and fabs(lmbda) < 1e273):
        return lgx
    else:
        return expm1(lmbda * lgx) / lmbda


cdef inline double inv_boxcox(double x, double lmbda) nogil:
    if lmbda == 0:
        return exp(x)
    else:
        return exp(log1p(lmbda * x) / lmbda)


cdef inline double inv_boxcox1p(double x, double lmbda) nogil:
    if lmbda == 0:
        return expm1(x)
    elif fabs(lmbda * x) < 1e-154:
        return x
    else:
        return expm1(log1p(lmbda * x) / lmbda)
