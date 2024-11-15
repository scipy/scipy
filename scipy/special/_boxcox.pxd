
from libc.math cimport log, log1p, expm1, exp, fabs, copysign


cdef inline double boxcox(double x, double lmbda) noexcept nogil:
    # if lmbda << 1 and log(x) < 1.0, the lmbda*log(x) product can lose
    # precision, furthermore, expm1(x) == x for x < eps.
    # For doubles, the range of log is -744.44 to +709.78, with eps being
    # the smallest value produced.  This range means that we will have
    # abs(lmbda)*log(x) < eps whenever abs(lmbda) <= eps/-log(min double)
    # which is ~2.98e-19.  
    if fabs(lmbda) < 1e-19:
        return log(x)
    elif lmbda * log(x) < 709.78:
        return expm1(lmbda * log(x)) / lmbda
    else:
        return copysign(1., lmbda) * exp(lmbda * log(x) - log(fabs(lmbda))) - 1 / lmbda


cdef inline double boxcox1p(double x, double lmbda) noexcept nogil:
    # The argument given above in boxcox applies here with the modification
    # that the smallest value produced by log1p is the minimum representable
    # value, rather than eps.  The second condition here prevents unflow
    # when log1p(x) is < eps.
    cdef double lgx = log1p(x)
    if fabs(lmbda) < 1e-19 or (fabs(lgx) < 1e-289 and fabs(lmbda) < 1e273):
        return lgx
    elif lmbda * lgx < 709.78:
        return expm1(lmbda * lgx) / lmbda
    else:
        return copysign(1., lmbda) * exp(lmbda * lgx - log(fabs(lmbda))) - 1 / lmbda


cdef inline double inv_boxcox(double x, double lmbda) noexcept nogil:
    if lmbda == 0:
        return exp(x)
    elif lmbda * x < 1.79e308:
        return exp(log1p(lmbda * x) / lmbda)
    else:
        return exp((log(copysign(1., lmbda) * (x + 1 / lmbda)) + log(fabs(lmbda))) / lmbda)


cdef inline double inv_boxcox1p(double x, double lmbda) noexcept nogil:
    if lmbda == 0:
        return expm1(x)
    elif fabs(lmbda * x) < 1e-154:
        return x
    elif lmbda * x < 1.79e308:
        return expm1(log1p(lmbda * x) / lmbda)
    else:
        return expm1((log(copysign(1., lmbda) * (x + 1 / lmbda)) + log(fabs(lmbda))) / lmbda)
