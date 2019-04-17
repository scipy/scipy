# -*- cython -*-
"""
Evaluate orthogonal polynomial values using recurrence relations.

References
----------

.. [AMS55] Abramowitz & Stegun, Section 22.5.

.. [MH] Mason & Handscombe, Chebyshev Polynomials, CRC Press (2003).

.. [LP] P. Levrie & R. Piessens, A note on the evaluation of orthogonal
        polynomials using recurrence relations, Internal Report TW74 (1985)
        Dept. of Computer Science, K.U. Leuven, Belgium
        https://lirias.kuleuven.be/handle/123456789/131600

"""
#
# Authors: Pauli Virtanen, Eric Moore
#

#------------------------------------------------------------------------------
# Direct evaluation of polynomials
#------------------------------------------------------------------------------
cimport cython
from libc.math cimport sqrt, exp, floor, fabs, log, sin, M_PI as pi

from numpy cimport npy_cdouble
from ._complexstuff cimport (
    nan, inf, number_t, npy_cdouble_from_double_complex,
    double_complex_from_npy_cdouble)

from . cimport sf_error
from ._cephes cimport Gamma, lgam, beta, lbeta
from ._cephes cimport hyp2f1 as hyp2f1_wrap

cdef extern from "specfun_wrappers.h":
    double hyp1f1_wrap(double a, double b, double x) nogil
    npy_cdouble chyp2f1_wrap( double a, double b, double c, npy_cdouble z) nogil
    npy_cdouble chyp1f1_wrap( double a, double b, npy_cdouble z) nogil

cdef extern from "c_misc/misc.h":
    double gammasgn(double x) nogil

# Fused type wrappers

cdef inline number_t hyp2f1(double a, double b, double c, number_t z) nogil:
    cdef npy_cdouble r
    if number_t is double:
        return hyp2f1_wrap(a, b, c, z)
    else:
        r = chyp2f1_wrap(a, b, c, npy_cdouble_from_double_complex(z))
        return double_complex_from_npy_cdouble(r)

cdef inline number_t hyp1f1(double a, double b, number_t z) nogil:
    cdef npy_cdouble r
    if number_t is double:
        return hyp1f1_wrap(a, b, z)
    else:
        r = chyp1f1_wrap(a, b, npy_cdouble_from_double_complex(z))
        return double_complex_from_npy_cdouble(r)

#-----------------------------------------------------------------------------
# Binomial coefficient
#-----------------------------------------------------------------------------

@cython.cdivision(True)
cdef inline double binom(double n, double k) nogil:
    cdef double kx, nx, num, den, dk, sgn
    cdef int i

    if n < 0:
        nx = floor(n)
        if n == nx:
            # undefined
            return nan

    kx = floor(k)
    if k == kx and (fabs(n) > 1e-8 or n == 0):
        # Integer case: use multiplication formula for less rounding error
        # for cases where the result is an integer.
        #
        # This cannot be used for small nonzero n due to loss of
        # precision.

        nx = floor(n)
        if nx == n and kx > nx/2 and nx > 0:
            # Reduce kx by symmetry
            kx = nx - kx

        if kx >= 0 and kx < 20:
            num = 1.0
            den = 1.0
            for i in range(1, 1 + <int>kx):
                num *= i + n - kx
                den *= i
                if fabs(num) > 1e50:
                    num /= den
                    den = 1.0
            return num/den

    # general case:
    if n >= 1e10*k and k > 0:
        # avoid under/overflows in intermediate results
        return exp(-lbeta(1 + n - k, 1 + k) - log(n + 1))
    elif k > 1e8*fabs(n):
        # avoid loss of precision
        num = Gamma(1 + n) / fabs(k) + Gamma(1 + n) * n / (2*k**2) # + ...
        num /= pi * fabs(k)**n
        if k > 0:
            kx = floor(k)
            if <int>kx == kx:
                dk = k - kx
                sgn = 1 if (<int>kx) % 2 == 0 else -1
            else:
                dk = k
                sgn = 1
            return num * sin((dk-n)*pi) * sgn
        else:
            kx = floor(k)
            if <int>kx == kx:
                return 0
            else:
                return num * sin(k*pi)
    else:
        return 1/(n + 1)/beta(1 + n - k, 1 + k)

#-----------------------------------------------------------------------------
# Jacobi
#-----------------------------------------------------------------------------

cdef inline number_t eval_jacobi(double n, double alpha, double beta, number_t x) nogil:
    cdef double a, b, c, d
    cdef number_t g

    d = binom(n+alpha, n)
    a = -n
    b = n + alpha + beta + 1
    c = alpha + 1
    g = 0.5*(1-x)
    return d * hyp2f1(a, b, c, g)

@cython.cdivision(True)
cdef inline double eval_jacobi_l(long n, double alpha, double beta, double x) nogil:
    cdef long kk
    cdef double p, d
    cdef double k, t

    if n < 0:
        return eval_jacobi(n, alpha, beta, x)
    elif n == 0:
        return 1.0
    elif n == 1:
        return 0.5*(2*(alpha+1)+(alpha+beta+2)*(x-1))
    else:
        d = (alpha+beta+2)*(x - 1) / (2*(alpha+1))
        p = d + 1
        for kk in range(n-1):
            k = kk+1.0
            t = 2*k+alpha+beta
            d = ((t*(t+1)*(t+2))*(x-1)*p + 2*k*(k+beta)*(t+2)*d) / (2*(k+alpha+1)*(k+alpha+beta+1)*t)
            p = d + p
        return binom(n+alpha, n)*p

#-----------------------------------------------------------------------------
# Shifted Jacobi
#-----------------------------------------------------------------------------

@cython.cdivision(True)
cdef inline number_t eval_sh_jacobi(double n, double p, double q, number_t x) nogil:
    return eval_jacobi(n, p-q, q-1, 2*x-1) / binom(2*n + p - 1, n)

@cython.cdivision(True)
cdef inline double eval_sh_jacobi_l(long n, double p, double q, double x) nogil:
    return eval_jacobi_l(n, p-q, q-1, 2*x-1) / binom(2*n + p - 1, n)

#-----------------------------------------------------------------------------
# Gegenbauer (Ultraspherical)
#-----------------------------------------------------------------------------

@cython.cdivision(True)
cdef inline number_t eval_gegenbauer(double n, double alpha, number_t x) nogil:
    cdef double a, b, c, d
    cdef number_t g

    d = Gamma(n+2*alpha)/Gamma(1+n)/Gamma(2*alpha)
    a = -n
    b = n + 2*alpha
    c = alpha + 0.5
    g = (1-x)/2.0
    return d * hyp2f1(a, b, c, g)

@cython.cdivision(True)
cdef inline double eval_gegenbauer_l(long n, double alpha, double x) nogil:
    cdef long kk
    cdef long a, b
    cdef double p, d
    cdef double k

    if n < 0:
        return 0.0
    elif n == 0:
        return 1.0
    elif n == 1:
        return 2*alpha*x
    elif alpha == 0.0:
        return eval_gegenbauer(n, alpha, x)
    elif fabs(x) < 1e-5:
        # Power series rather than recurrence due to loss of precision
        # http://functions.wolfram.com/Polynomials/GegenbauerC3/02/
        a = n//2

        d = 1 if a % 2 == 0 else -1
        d /= beta(alpha, 1 + a)
        if n == 2*a:
            d /= (a + alpha)
        else:
            d *= 2*x

        p = 0
        for kk in range(a+1):
            p += d
            d *= -4*x**2 * (a - kk) * (-a + alpha + kk + n) / (
                (n + 1 - 2*a + 2*kk) * (n + 2 - 2*a + 2*kk))
            if fabs(d) == 1e-20*fabs(p):
                # converged
                break
        return p
    else:
        d = x - 1
        p = x
        for kk in range(n-1):
            k = kk+1.0
            d = (2*(k+alpha)/(k+2*alpha))*(x-1)*p + (k/(k+2*alpha)) * d
            p = d + p

        if fabs(alpha/n) < 1e-8:
            # avoid loss of precision
            return 2*alpha/n * p
        else:
            return binom(n+2*alpha-1, n)*p

#-----------------------------------------------------------------------------
# Chebyshev 1st kind (T)
#-----------------------------------------------------------------------------

cdef inline number_t eval_chebyt(double n, number_t x) nogil:
    cdef double a, b, c, d
    cdef number_t g

    d = 1.0
    a = -n
    b = n
    c = 0.5
    g = 0.5*(1-x)
    return hyp2f1(a, b, c, g)

cdef inline double eval_chebyt_l(long k, double x) nogil:
    # Use Chebyshev T recurrence directly, see [MH]
    cdef long m
    cdef double b2, b1, b0

    if k < 0:
        # symmetry
        k = -k

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

cdef inline number_t eval_chebyu(double n, number_t x) nogil:
    cdef double a, b, c, d
    cdef number_t g

    d = n+1
    a = -n
    b = n+2
    c = 1.5
    g = 0.5*(1-x)
    return d*hyp2f1(a, b, c, g)

cdef inline double eval_chebyu_l(long k, double x) nogil:
    cdef long m
    cdef int sign
    cdef double b2, b1, b0

    if k == -1:
        return 0
    elif k < -1:
        # symmetry
        k = -k - 2
        sign = -1
    else:
        sign = 1

    b2 = 0
    b1 = -1
    b0 = 0
    x = 2*x
    for m in range(k+1):
        b2 = b1
        b1 = b0
        b0 = x*b1 - b2
    return b0 * sign

#-----------------------------------------------------------------------------
# Chebyshev S
#-----------------------------------------------------------------------------

cdef inline number_t eval_chebys(double n, number_t x) nogil:
    return eval_chebyu(n, 0.5*x)

cdef inline double eval_chebys_l(long n, double x) nogil:
    return eval_chebyu_l(n, 0.5*x)

#-----------------------------------------------------------------------------
# Chebyshev C
#-----------------------------------------------------------------------------

cdef inline number_t eval_chebyc(double n, number_t x) nogil:
    return 2*eval_chebyt(n, 0.5*x)

cdef inline double eval_chebyc_l(long n, double x) nogil:
    return 2*eval_chebyt_l(n, 0.5*x)

#-----------------------------------------------------------------------------
# Chebyshev 1st kind shifted
#-----------------------------------------------------------------------------

cdef inline number_t eval_sh_chebyt(double n, number_t x) nogil:
    return eval_chebyt(n, 2*x-1)

cdef inline double eval_sh_chebyt_l(long n, double x) nogil:
    return eval_chebyt_l(n, 2*x-1)

#-----------------------------------------------------------------------------
# Chebyshev 2st kind shifted
#-----------------------------------------------------------------------------

cdef inline number_t eval_sh_chebyu(double n, number_t x) nogil:
    return eval_chebyu(n, 2*x-1)

cdef inline double eval_sh_chebyu_l(long n, double x) nogil:
    return eval_chebyu_l(n, 2*x-1)

#-----------------------------------------------------------------------------
# Legendre
#-----------------------------------------------------------------------------

cdef inline number_t eval_legendre(double n, number_t x) nogil:
    cdef double a, b, c, d
    cdef number_t g

    d = 1
    a = -n
    b = n+1
    c = 1
    g = 0.5*(1-x)
    return d*hyp2f1(a, b, c, g)

@cython.cdivision(True)
cdef inline double eval_legendre_l(long n, double x) nogil:
    cdef long kk, a
    cdef double p, d
    cdef double k

    if n < 0:
        # symmetry
        n = -n - 1

    if n == 0:
        return 1.0
    elif n == 1:
        return x
    elif fabs(x) < 1e-5:
        # Power series rather than recurrence due to loss of precision
        # http://functions.wolfram.com/Polynomials/LegendreP/02/
        a = n//2

        d = 1 if a % 2 == 0 else -1
        if n == 2*a:
            d *= -2 / beta(a + 1, -0.5)
        else:
            d *= 2 * x / beta(a + 1, 0.5)

        p = 0
        for kk in range(a+1):
            p += d
            d *= -2 * x**2 * (a - kk) * (2*n + 1 - 2*a + 2*kk) / (
                (n + 1 - 2*a + 2*kk) * (n + 2 - 2*a + 2*kk))
            if fabs(d) == 1e-20*fabs(p):
                # converged
                break
        return p
    else:
        d = x - 1
        p = x
        for kk in range(n-1):
            k = kk+1.0
            d = ((2*k+1)/(k+1))*(x-1)*p + (k/(k+1)) * d
            p = d + p
        return p

#-----------------------------------------------------------------------------
# Legendre Shifted
#-----------------------------------------------------------------------------

cdef inline number_t eval_sh_legendre(double n, number_t x) nogil:
    return eval_legendre(n, 2*x-1)

cdef inline double eval_sh_legendre_l(long n, double x) nogil:
    return eval_legendre_l(n, 2*x-1)

#-----------------------------------------------------------------------------
# Generalized Laguerre
#-----------------------------------------------------------------------------

cdef inline number_t eval_genlaguerre(double n, double alpha, number_t x) nogil:
    cdef double a, b, d
    cdef number_t g

    if alpha <= -1:
        sf_error.error("eval_genlaguerre", sf_error.DOMAIN,
                       "polynomial defined only for alpha > -1")
        return nan

    d = binom(n+alpha, n)
    a = -n
    b = alpha + 1
    g = x
    return d * hyp1f1(a, b, g)

@cython.cdivision(True)
cdef inline double eval_genlaguerre_l(long n, double alpha, double x) nogil:
    cdef long kk
    cdef double p, d
    cdef double k

    if alpha <= -1:
        sf_error.error("eval_genlaguerre", sf_error.DOMAIN,
                       "polynomial defined only for alpha > -1")
        return nan

    if n < 0:
        return 0.0
    elif n == 0:
        return 1.0
    elif n == 1:
        return -x+alpha+1
    else:
        d = -x/(alpha+1)
        p = d + 1
        for kk in range(n-1):
            k = kk+1.0
            d = -x/(k+alpha+1)*p + (k/(k+alpha+1)) * d
            p = d + p
        return binom(n+alpha, n)*p

#-----------------------------------------------------------------------------
# Laguerre
#-----------------------------------------------------------------------------

cdef inline number_t eval_laguerre(double n, number_t x) nogil:
    return eval_genlaguerre(n, 0., x)

cdef inline double eval_laguerre_l(long n, double x) nogil:
    return eval_genlaguerre_l(n, 0., x)

#-----------------------------------------------------------------------------
# Hermite (statistician's)
#-----------------------------------------------------------------------------

cdef inline double eval_hermitenorm(long n, double x) nogil:
    cdef long k
    cdef double y1, y2, y3

    if n < 0:
        return 0.0
    elif n == 0:
        return 1.0
    elif n == 1:
        return x
    else:
        y3 = 0.0
        y2 = 1.0
        for k in range(n, 1, -1):
            y1 = x*y2 - k*y3
            y3 = y2
            y2 = y1
        return x*y2 - y3

#-----------------------------------------------------------------------------
# Hermite (physicist's)
#-----------------------------------------------------------------------------

@cython.cdivision(True)
cdef inline double eval_hermite(long n, double x) nogil:
    return eval_hermitenorm(n, sqrt(2)*x) * 2**(n/2.0)
