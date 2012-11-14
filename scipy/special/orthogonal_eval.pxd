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
cimport cython
from libc.math cimport sqrt, exp

cdef extern from "numpy/npy_math.h":
    ctypedef struct npy_cdouble:
        double real
        double imag

cdef extern from "cephes.h":
    double Gamma(double x) nogil
    double lgam(double x) nogil
    double hyp2f1 (double a, double b, double c, double x) nogil 

cdef extern from "specfun_wrappers.h":
    double hyp1f1_wrap(double a, double b, double x) nogil
    npy_cdouble chyp2f1_wrap( double a, double b, double c, npy_cdouble z) nogil 
    npy_cdouble chyp1f1_wrap( double a, double b, npy_cdouble z) nogil
  

cdef inline double binom(double n, double k):
    return exp(lgam(n+1) - lgam(k+1) - lgam(1+n-k))

#-----------------------------------------------------------------------------
# Jacobi
#-----------------------------------------------------------------------------

cdef inline double eval_jacobi_dddd(double n, double alpha, double beta, double x):
    cdef double a, b, c, d, g
    
    d = binom(n+alpha, n)
    a = -n
    b = n + alpha + beta + 1
    c = alpha + 1
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

cdef inline double complex eval_jacobi_dddD(double n, double alpha, double beta, double complex x):
    cdef double a, b, c, d 
    cdef double complex g
    cdef npy_cdouble r
    
    d = binom(n+alpha, n)
    a = -n
    b = n + alpha + beta + 1
    c = alpha + 1
    g = (1-x)/2.0
    r = chyp2f1_wrap(a, b, c, (<npy_cdouble*>&g)[0])
    r.real *= d
    r.imag *= d
    return (<double complex*>&r)[0]

#-----------------------------------------------------------------------------
# Shifted Jacobi
#-----------------------------------------------------------------------------

cdef inline double eval_sh_jacobi_dddd(double n, double p, double q, double x):
    cdef double factor

    factor = exp(lgam(1+n) + lgam(n+p) - lgam(2*n+p))
    return factor * eval_jacobi_dddd(n, p-q, q-1, 2*x-1)

cdef inline double complex eval_sh_jacobi_dddD(double n, double p, double q, double complex x):
    cdef double factor

    factor = exp(lgam(1+n) + lgam(n+p) - lgam(2*n+p))
    return factor * eval_jacobi_dddD(n, p-q, q-1, 2*x-1) 

#-----------------------------------------------------------------------------
# Gegenbauer (Ultraspherical)
#-----------------------------------------------------------------------------

cdef inline double eval_gegenbauer_ddd(double n, double alpha, double x):
    cdef double a, b, c, d, g

    d = Gamma(n+2*alpha)/Gamma(1+n)/Gamma(2*alpha)
    a = -n
    b = n + 2*alpha
    c = alpha + 0.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

cdef inline double complex eval_gegenbauer_ddD(double n, double alpha, double complex x):
    cdef double a, b, c, d
    cdef double complex g
    cdef npy_cdouble r

    d = Gamma(n+2*alpha)/Gamma(1+n)/Gamma(2*alpha)
    a = -n
    b = n + 2*alpha
    c = alpha + 0.5
    g = (1-x)/2.0
    r = chyp2f1_wrap(a, b, c, (<npy_cdouble*>&g)[0])
    r.real *= d
    r.imag *= d
    return (<double complex*>&r)[0]

#-----------------------------------------------------------------------------
# Chebyshev 1st kind (T)
#-----------------------------------------------------------------------------

cdef inline double eval_chebyt_dd(double n, double x) nogil:
    cdef double a, b, c, d, g

    d = 1.0
    a = -n
    b = n 
    c = 0.5
    g = (1-x)/2.0
    return hyp2f1(a,b,c,g)

cdef inline double complex eval_chebyt_dD(double n, double complex x):
    cdef double a, b, c, d
    cdef double complex g
    cdef npy_cdouble r

    d = 1.0
    a = -n
    b = n
    c = 0.5
    g = (1-x)/2.0
    r = chyp2f1_wrap(a, b, c, (<npy_cdouble*>&g)[0])
    return (<double complex*>&r)[0]

cdef inline double eval_chebyt_ld(long k, double x) nogil:
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

#-----------------------------------------------------------------------------
# Chebyshev 2st kind (U)
#-----------------------------------------------------------------------------

cdef inline double eval_chebyu_dd(double n, double x):
    cdef double a, b, c, d, g

    d = n+1
    a = -n
    b = n+2
    c = 1.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

cdef inline double complex eval_chebyu_dD(double n, double complex x):
    cdef double a, b, c, d
    cdef double complex g
    cdef npy_cdouble r

    d = n+1
    a = -n
    b = n+2
    c = 1.5
    g = (1-x)/2.0
    r = chyp2f1_wrap(a, b, c, (<npy_cdouble*>&g)[0])
    r.real *= d
    r.imag *= d
    return (<double complex*>&r)[0]

#-----------------------------------------------------------------------------
# Chebyshev S
#-----------------------------------------------------------------------------

cdef inline double eval_chebys_dd(double n, double x):
    return eval_chebyu_dd(n, x/2.0)

cdef inline double complex eval_chebys_dD(double n, double complex x):
    return eval_chebyu_dD(n, x/2.0)

#-----------------------------------------------------------------------------
# Chebyshev C
#-----------------------------------------------------------------------------

cdef inline double eval_chebyc_dd(double n, double x):
    return 2*eval_chebyt_dd(n, x/2.0)

cdef inline double complex eval_chebyc_dD(double n, double complex x):
    return 2*eval_chebyt_dD(n, x/2.0)

cdef inline double eval_chebyc_ld(long n, double x):
    return 2*eval_chebyt_ld(n, x/2.0)

#-----------------------------------------------------------------------------
# Chebyshev 1st kind shifted
#-----------------------------------------------------------------------------

cdef inline double eval_sh_chebyt_dd(double n, double x):
    return eval_chebyt_dd(n, 2*x-1)

cdef inline double complex eval_sh_chebyt_dD(double n, double complex x):
    return eval_chebyt_dD(n, 2*x-1)

cdef inline double eval_sh_chebyt_ld(long n, double x):
    return eval_chebyt_ld(n, 2*x-1)

#-----------------------------------------------------------------------------
# Chebyshev 2st kind shifted
#-----------------------------------------------------------------------------

cdef inline double eval_sh_chebyu_dd(double n, double x):
    return eval_chebyu_dd(n, 2*x-1)

cdef inline double complex eval_sh_chebyu_dD(double n, double complex x):
    return eval_chebyu_dD(n, 2*x-1)

#-----------------------------------------------------------------------------
# Legendre
#-----------------------------------------------------------------------------

cdef inline double eval_legendre_dd(double n, double x) nogil:
    cdef double  a, b, c, d, g

    d = 1
    a = -n
    b = n+1
    c = 1
 #   g = (1-x)/2.0
    g = 0.5 * (1-x)
    return hyp2f1(a, b, c, g) * d

cdef inline double complex eval_legendre_dD(double n, double complex x):
    cdef double a, b, c, d
    cdef double complex g
    cdef npy_cdouble r

    d = 1
    a = -n
    b = n+1
    c = 1
    g = (1-x)/2.0
    r = chyp2f1_wrap(a, b, c, (<npy_cdouble*>&g)[0])
    r.real *= d
    r.imag *= d
    return (<double complex*>&r)[0]

#-----------------------------------------------------------------------------
# Legendre Shifted
#-----------------------------------------------------------------------------

cdef inline double eval_sh_legendre_dd(double n, double x):
    return eval_legendre_dd(n, 2*x-1)

cdef inline double complex eval_sh_legendre_dD(double n, double complex x):
    return eval_legendre_dD(n, 2*x-1)

#-----------------------------------------------------------------------------
# Generalized Laguerre
#-----------------------------------------------------------------------------

cdef inline double eval_genlaguerre_ddd(double n, double alpha, double x):
    cdef double a, b, d, g

    d = binom(n+alpha, n)
    a = -n
    b = alpha + 1
    g = x
    return hyp1f1_wrap(a, b, g) * d

cdef inline double complex eval_genlaguerre_ddD(double n, double alpha, double complex x):
    cdef double a, b, d
    cdef double complex g
    cdef npy_cdouble r

    d = binom(n+alpha, n)
    a = -n
    b = alpha + 1
    g = x
    r = chyp1f1_wrap(a, b, (<npy_cdouble*>&g)[0])
    r.real *= d
    r.imag *= d
    return (<double complex*>&r)[0]

#-----------------------------------------------------------------------------
# Laguerre
#-----------------------------------------------------------------------------

cdef inline double eval_laguerre_dd(double n, double x):
    return eval_genlaguerre_ddd(n, 0., x)

cdef inline double complex eval_laguerre_dD(double n, double complex x):
    return eval_genlaguerre_ddD(n, 0., x)

#-----------------------------------------------------------------------------
# Hermite (physicist's)
#-----------------------------------------------------------------------------

cdef inline double eval_hermite(long n, double x):
    cdef long m

    if n % 2 == 0:
        m = n/2
        return ((-1)**m * 2**(2*m) * Gamma(1+m)
                 * eval_genlaguerre_ddd(m, -0.5, x**2))
    else:
        m = (n-1)/2
        return ((-1)**m * 2**(2*m+1) * Gamma(1+m)
                  * x * eval_genlaguerre_ddd(m, 0.5, x**2))

#-----------------------------------------------------------------------------
# Hermite (statistician's)
#-----------------------------------------------------------------------------

cdef inline double eval_hermitenorm(long n, double x):
    return eval_hermite(n, x/sqrt(2)) * 2**(-n/2.0)

