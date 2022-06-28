# Implementation of sin/cos/sinh/cosh integrals for complex arguments
#
# Sources
# [1] Fredrik Johansson and others. mpmath: a Python library for
#     arbitrary-precision floating-point arithmetic (version 0.19),
#     December 2013. http://mpmath.org/.
# [2] NIST, "Digital Library of Mathematical Functions",
#     https://dlmf.nist.gov/
import cython
cimport numpy as np

from . cimport sf_error

cdef extern from "numpy/npy_math.h":
    double NPY_PI
    double NPY_PI_2
    double NPY_EULER

from ._complexstuff cimport (
    npy_cdouble_from_double_complex, double_complex_from_npy_cdouble,
    inf, nan, zabs, zlog, zpack)

cdef extern from "specfun_wrappers.h":
    np.npy_cdouble cexpi_wrap(np.npy_cdouble) nogil

DEF MAXITER = 100
DEF TOL = 2.220446092504131e-16
    

cdef inline double complex zexpi(double complex z) nogil:
    cdef np.npy_cdouble r
    r = cexpi_wrap(npy_cdouble_from_double_complex(z))
    return double_complex_from_npy_cdouble(r)


@cython.cdivision(True)
cdef inline void power_series(int sgn, double complex z,
                             double complex *s, double complex *c) nogil:
    """DLMF 6.6.5 and 6.6.6. If sgn = -1 computes si/ci, and if sgn = 1
    computes shi/chi.

    """
    cdef:
        int n
        double complex fac, term1, term2
        
    fac = z
    s[0] = fac
    c[0] = 0        
    for n in range(1, MAXITER):
        fac *= sgn*z/(2*n)
        term2 = fac/(2*n)
        c[0] += term2
        fac *= z/(2*n + 1)
        term1 = fac/(2*n + 1)
        s[0] += term1
        if zabs(term1) < TOL*zabs(s[0]) and zabs(term2) < TOL*zabs(c[0]):
            break

    
cdef inline int csici(double complex z,
                      double complex *si, double complex *ci) nogil:
    """Compute sin/cos integrals at complex arguments. The algorithm
    largely follows that of [1].

    """
    cdef double complex jz, term1, term2

    if z == inf:
        si[0] = NPY_PI_2
        ci[0] = 0
        return 0
    elif z == -inf:
        si[0] = -NPY_PI_2
        ci[0] = 1j*NPY_PI
        return 0
    elif zabs(z) < 0.8:
        # Use the series to avoid cancellation in si
        power_series(-1, z, si, ci)
        if z == 0:
            sf_error.error("sici", sf_error.DOMAIN, NULL)
            ci[0] = zpack(-inf, nan)
        else:
            ci[0] += NPY_EULER + zlog(z)
        return 0
    
    # DLMF 6.5.5/6.5.6 plus DLMF 6.4.4/6.4.6/6.4.7
    jz = 1j*z
    term1 = zexpi(jz)
    term2 = zexpi(-jz)
    si[0] = -0.5j*(term1 - term2)
    ci[0] = 0.5*(term1 + term2)
    if z.real == 0:
        if z.imag > 0:
            ci[0] += 1j*NPY_PI_2
        elif z.imag < 0:
            ci[0] -= 1j*NPY_PI_2
    elif z.real > 0:
        si[0] -= NPY_PI_2
    else:
        si[0] += NPY_PI_2
        if z.imag >= 0:
            ci[0] += 1j*NPY_PI
        else:
            ci[0] -= 1j*NPY_PI

    return 0


cdef inline int cshichi(double complex z,
                        double complex *shi, double complex *chi) nogil:
    """Compute sinh/cosh integrals at complex arguments. The algorithm
    largely follows that of [1].

    """
    cdef double complex term1, term2

    if z == inf:
        shi[0] = inf
        chi[0] = inf
        return 0
    elif z == -inf:
        shi[0] = -inf
        chi[0] = inf
        return 0
    elif zabs(z) < 0.8:
        # Use the series to avoid cancellation in shi
        power_series(1, z, shi, chi)
        if z == 0:
            sf_error.error("shichi", sf_error.DOMAIN, NULL)
            chi[0] = zpack(-inf, nan)
        else:
            chi[0] += NPY_EULER + zlog(z)
        return 0

    term1 = zexpi(z)
    term2 = zexpi(-z)
    shi[0] = 0.5*(term1 - term2)
    chi[0] = 0.5*(term1 + term2)
    if z.imag > 0:
        shi[0] -= 0.5j*NPY_PI
        chi[0] += 0.5j*NPY_PI
    elif z.imag < 0:
        shi[0] += 0.5j*NPY_PI
        chi[0] -= 0.5j*NPY_PI
    elif z.real < 0:
        chi[0] += 1j*NPY_PI

    return 0
