# Compute the exponential integral.

from ._floatingstuff cimport *


cdef inline floating expi(floating x) nogil:
    """Compute the exponential integral Ei(x)."""
    cdef int k
    cdef floating ei, one, r

    if x == <floating>0:
        ei = -infinity(x)
    elif x < <floating>0:
        ei = -exp1(-x)
    elif fabs(x) <= <floating>40:
        # Power series around x = 0
        ei = <floating>1
        r = <floating>1
        for k in range(1, 101):
            r *= k*x/(k + <floating>1)**2
            ei += r
            if fabs(r) <= eps(x)*fabs(ei):
                break
        ei = euler(x) + log(x) + x*ei
    else:
        # Asymptotic expansion
        ei = <floating>1
        r = <floating>1
        for k in range(1, 21):
            r *= k/x
            ei += r
        ei = exp(x)/x*ei
    return ei


cdef inline floating exp1(floating x) nogil:
    """Compute the exponential integral E1(x)."""
    cdef int k, m
    cdef floating e1, r, t, t0

    if x == <floating>0:
        e1 = infinity(x)
    if x < <floating>0:
        return nan(x)
    elif x <= <floating>1:
        e1 = <floating>1
        r = <floating>1
        for k in range(1,26):
            r *= -k*x/(k + <floating>1)**2
            e1 += r
            if fabs(r) <= eps(x)*fabs(e1):
                break
        e1 = -euler(x) - log(x) + x*e1
    else:
        m = 20 + <int>(<floating>80/x)
        t0 = <floating>0
        # Work around https://github.com/cython/cython/issues/532
        k = m
        while k > 0:
            t0 = k/(<floating>1 + k/(x+t0))
            k -= 1
        t = (<floating>1)/(x + t0)
        e1 = exp(-x)*t
    return e1
