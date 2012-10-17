"""
Evaluate orthogonal polynomial values using recurrence relations.

References
----------

.. [AMS55] Abramowitz & Stegun, Section 22.5.

.. [MH] Mason & Handscombe, Chebyshev Polynomials, CRC Press (2003).

"""
#
# Copyright (C) 2009 Pauli Virtanen
# Distributed under the same license as Scipy.
#

#------------------------------------------------------------------------------
# Direct evaluation of polynomials
#------------------------------------------------------------------------------

from libc.math cimport sqrt, exp

cdef extern from "cephes.h":
    double Gamma(double x)
    double lgam(double x)
    double hyp2f1 (double a, double b, double c, double x) 

cdef extern from "amos_wrappers.h":
    double hyp1f1_wrap(double a, double b, double x)

cdef inline double binom(double n, double k):
    return exp(lgam(n+1) - lgam(k+1) - lgam(1+n-k))

cdef inline double eval_jacobi(double n, double alpha, double beta, double x):
    cdef double a, b, c, d, g
    
    d = binom(n+alpha, n)
    a = -n
    b = n + alpha + beta + 1
    c = alpha + 1
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

cdef inline double eval_sh_jacobi(double n, double p, double q, double x):
    cdef double factor

    factor = exp(lgam(1+n) + lgam(n+p) - lgam(2*n+p))
    return factor * eval_jacobi(n, p-q, q-1, 2*x-1)

cdef inline double eval_gegenbauer(double n, double alpha, double x):
    cdef double a, b, c, d, g

    d = Gamma(n+2*alpha)/Gamma(1+n)/Gamma(2*alpha)
    a = -n
    b = n + 2*alpha
    c = alpha + 0.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

cdef inline double eval_chebyt(long k, double x) nogil:
    # Use Chebyshev T recurrence directly, see [MH]
    cdef long m
    cdef double b2, b1, b0

    b2 = 0
    b1 = -1
    b0 = 0
    x = 2*x
    for m in range(k+1):
        b2 = b1
        b1 = b0
        b0 = x*b1 - b2
    return (b0 - b2)/2.0

cdef inline double eval_chebyu(double n, double x):
    cdef double a, b, c, d, g

    d = n+1
    a = -n
    b = n+2
    c = 1.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

cdef inline double eval_chebys(double n, double x):
    return eval_chebyu(n, x/2.0)

cdef inline double eval_chebyc(long n, double x):
    return 2*eval_chebyt(n, x/2.0)

cdef inline double eval_sh_chebyt(long n, double x):
    return eval_chebyt(n, 2*x-1)

cdef inline double eval_sh_chebyu(double n, double x):
    return eval_chebyu(n, 2*x-1)

cdef inline double eval_legendre(double n, double x):
    cdef double  a, b, c, d, g

    d = 1
    a = -n
    b = n+1
    c = 1
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

cdef inline double eval_sh_legendre(double n, double x):
    return eval_legendre(n, 2*x-1)

cdef inline double eval_genlaguerre(double n, double alpha, double x):
    cdef double a, b, d, g

    d = binom(n+alpha, n)
    a = -n
    b = alpha + 1
    g = x
    return hyp1f1_wrap(a, b, g) * d

cdef inline double eval_laguerre(double n, double x):
    return eval_genlaguerre(n, 0., x)

cdef inline double eval_hermite(long n, double x):
    cdef long m

    if n % 2 == 0:
        m = n/2
        return ((-1)**m * 2**(2*m) * Gamma(1+m)
                 * eval_genlaguerre(m, -0.5, x**2))
    else:
        m = (n-1)/2
        return ((-1)**m * 2**(2*m+1) * Gamma(1+m)
                  * x * eval_genlaguerre(m, 0.5, x**2))

cdef inline double eval_hermitenorm(long n, double x):
    return eval_hermite(n, x/sqrt(2)) * 2**(-n/2.0)

