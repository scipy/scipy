"""Evaluate polynomials.

All of the coefficients are stored in reverse order, i.e. if the
polynomial is

    u_n x^n + u_{n - 1} x^{n - 1} + ... + u_0,

then coeffs[0] = u_n, coeffs[1] = u_{n - 1}, ..., coeffs[n] = u_0.

References
----------
[1] Knuth, "The Art of Computer Programming, Volume II"

"""
from ._complexstuff cimport zabs

cdef extern from "_c99compat.h":
    double sc_fma(double x, double y, double z) nogil


cdef inline double complex cevalpoly(double *coeffs, int degree,
                                     double complex z) nogil:
    """Evaluate a polynomial with real coefficients at a complex point.

    Uses equation (3) in section 4.6.4 of [1]. Note that it is more
    efficient than Horner's method.

    """
    cdef:
        int j
        double a = coeffs[0]
        double b = coeffs[1]
        double r = 2*z.real
        double s = z.real*z.real + z.imag*z.imag
        double tmp

    for j in range(2, degree + 1):
        tmp = b
        b = sc_fma(-s, a, coeffs[j])
        a = sc_fma(r, a, tmp)
    return z*a + b
