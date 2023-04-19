from libc.math cimport log, fabs, expm1, log1p, isnan, NAN, INFINITY

cdef inline double entr(double x) noexcept nogil:
    if isnan(x):
        return x
    elif x > 0:
        return -x * log(x)
    elif x == 0:
        return 0
    else:
        return -INFINITY

cdef inline double kl_div(double x, double y) noexcept nogil:
    if isnan(x) or isnan(y):
        return NAN
    elif x > 0 and y > 0:
        return x * log(x / y) - x + y
    elif x == 0 and y >= 0:
        return y
    else:
        return INFINITY

cdef inline double rel_entr(double x, double y) noexcept nogil:
    if isnan(x) or isnan(y):
        return NAN
    elif x > 0 and y > 0:
        return x * log(x / y)
    elif x == 0 and y >= 0:
        return 0
    else:
        return INFINITY

cdef inline double huber(double delta, double r) noexcept nogil:
    if delta < 0:
        return INFINITY
    elif fabs(r) <= delta:
        return 0.5 * r * r;
    else:
        return delta * (fabs(r) - 0.5 * delta);

cdef inline double pseudo_huber(double delta, double r) noexcept nogil:
    cdef double u, v
    if delta < 0:
        return INFINITY
    elif delta == 0 or r == 0:
        return 0
    else:
        u = delta
        v = r / delta
        # The formula is u*u*(sqrt(1 + v*v) - 1), but to maintain
        # precision with small v, we use
        #   sqrt(1 + v*v) - 1  =  exp(0.5*log(1 + v*v)) - 1
        #                      =  expm1(0.5*log1p(v*v))
        return u*u*expm1(0.5*log1p(v*v))
