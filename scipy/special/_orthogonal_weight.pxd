"""Weight functions for orthogonal polynomials.

References
----------

.. [AMS55] Abramowitz & Stegun, Section 22.5.

"""
cimport cython
from libc.math cimport sqrt, exp

from numpy.math cimport NAN

cimport sf_error


@cython.cdivision(True)
cdef inline double weight_jacobi(double alpha, double beta, double x) nogil:
    if x < -1 or x > 1 or alpha <= -1 or beta <= -1:
        sf_error.error("weight_jacobi", sf_error.DOMAIN, NULL)
        return NAN

    return (1 - x)**alpha*(1 + x)**beta


@cython.cdivision(True)
cdef inline double weight_sh_jacobi(double p, double q, double x) nogil:
    if x < 0 or x > 1 or p - q <= -1 or q <= 0:
        sf_error.error("weight_sh_jacobi", sf_error.DOMAIN, NULL)
        return NAN

    return (1 - x)**(p - q)*x**(q - 1)


@cython.cdivision(True)
cdef inline double weight_gegenbauer(double alpha, double x) nogil:
    if x < -1 or x > 1 or alpha <= -0.5:
        sf_error.error("weight_gegenbauer", sf_error.DOMAIN, NULL)
        return NAN

    return (1 - x**2)**(alpha - 0.5)


@cython.cdivision(True)
cdef inline double weight_chebyt(double x) nogil:
    if x < -1 or x > 1:
        sf_error.error("weight_chebyt", sf_error.DOMAIN, NULL)
        return NAN

    return 1.0/sqrt(1 - x**2)


@cython.cdivision(True)
cdef inline double weight_chebyu(double x) nogil:
    if x < -1 or x > 1:
        sf_error.error("weight_chebyu", sf_error.DOMAIN, NULL)

    return sqrt(1 - x**2)


@cython.cdivision(True)
cdef inline double weight_chebyc(double x) nogil:
    if x < -2 or x > 2:
        sf_error.error("weight_chebyc", sf_error.DOMAIN, NULL)

    return 1.0/sqrt(1 - 0.25*x**2)


@cython.cdivision(True)
cdef inline double weight_chebys(double x) nogil:
    if x < -2 or x > 2:
        sf_error.error("weight_chebys", sf_error.DOMAIN, NULL)

    return sqrt(1 - 0.25*x**2)


@cython.cdivision(True)
cdef inline double weight_sh_chebyt(double x) nogil:
    if x < 0 or x > 1:
        sf_error.error("weight_sh_chebyt", sf_error.DOMAIN, NULL)

    return 1.0/sqrt(x - x**2)


@cython.cdivision(True)
cdef inline double weight_sh_chebyu(double x) nogil:
    if x < 0 or x > 1:
        sf_error.error("weight_sh_chebyu", sf_error.DOMAIN, NULL)

    return sqrt(x - x**2)


@cython.cdivision(True)
cdef inline double weight_legendre(double x) nogil:
    if x < -1 or x > 1:
        sf_error.error("weight_legendre", sf_error.DOMAIN, NULL)

    return 1.0


@cython.cdivision(True)
cdef inline double weight_sh_legendre(double x) nogil:
    if x < 0 or x > 1:
        sf_error.error("weight_sh_legendre", sf_error.DOMAIN, NULL)

    return 1.0


@cython.cdivision(True)
cdef inline double weight_genlaguerre(double alpha, double x) nogil:
    if x < 0 or alpha <= -1:
        sf_error.error("weight_genlaguerre", sf_error.DOMAIN, NULL)

    return exp(-x)*x**alpha


@cython.cdivision(True)
cdef inline double weight_laguerre(double x) nogil:
    if x < 0:
        sf_error.error("weight_laguerre", sf_error.DOMAIN, NULL)

    return exp(-x)


@cython.cdivision(True)
cdef inline double weight_hermite(double x) nogil:
    return exp(-x**2)


@cython.cdivision(True)
cdef inline double weight_hermitenorm(double x) nogil:
    return exp(-x**2/2)
