cimport cython
from numpy cimport npy_cdouble
from libc.math cimport fabs, floor, exp, M_LN2, M_PI, pow

from . cimport sf_error
from ._cephes cimport Gamma, gammasgn, lgam
from ._complexstuff cimport (
    double_complex_from_npy_cdouble, npy_cdouble_from_double_complex,
    zabs, zisinf, zisnan, zpack, zpow
)


cdef extern from "numpy/npy_math.h":
    double NPY_NAN
    double NPY_INFINITY


cdef extern from 'specfun_wrappers.h':
    npy_cdouble chyp2f1_wrap(
        double,
        double,
        double,
        npy_cdouble
    ) nogil


DEF EPS = 2.220446049250313e-16
DEF SQRT_PI = 1.7724538509055159  # sqrt(M_PI)
DEF LOG_PI_2 = 0.5723649429247001  # log(M_PI) / 2


@cython.cdivision(True)
cdef inline double complex hyp2f1_complex(
        double a, double b, double c, double complex z
) nogil:
    cdef:
        int max_iter
        double modulus_z
        double complex result
        bint condition1, condition2
    modulus_z = zabs(z)
    # Special Cases
    # -------------------------------------------------------------------------
    # Diverges when c is a negative integer
    if c == floor(c) and c < 0:
        return NPY_INFINITY + 0.0j
    # Diverges as real(z) -> 1 when c < a + b.
    if fabs(1 - z.real) < EPS and z.imag == 0 and c - a - b < 0:
        return NPY_INFINITY + 0.0j
    # Equals 1 at z = 0. Constant function 1 when a = 0 or b = 0.
    if modulus_z == 0 or a == 0 or b == 0:
        return 1.0 + 0.0j
    # Gauss's Summation Theorem for z = 1; c - a - b > 0 (DLMF 15.4.20).
    if z == 1.0 and c - a - b > 0:
        result = Gamma(c) * Gamma(c - a - b)
        result /= Gamma(c - a) * Gamma(c - b)
        # If there has been an overflow, try again with logs
        if zisnan(result) or zisinf(result) or result == 0.0:
            result = exp(
                lgam(c) - lgam(c - a) +
                lgam(c - a - b) - lgam(c - b)
            )
            result *= gammasgn(c) * gammasgn(c - a - b)
            result *= gammasgn(c - a) * gammasgn(c - b)
        return result
    # Kummer's Theorem for z = -1; c = 1 + a - b (DLMF 15.4.26).
    # hyp2f1(a, b, 1 + a - b, -1) =
    # gamma(1 + a - b) * gamma(1 + 0.5 * a) /
    # (gamma(1 + a) * gamma(1 + 0.5*a - b)
    if zabs(z + 1) < EPS and fabs(1 + a - b - c) < EPS:
        # The computation below has been simplified through
        # Legendre duplication for the Gamma function (DLMF 5.5.5).
        # gamma(1 + 0.5*a)*gamma(0.5 + 0.5*a) = sqrt(pi)*2**(-a)*gamma(1 + a)
        result = SQRT_PI * pow(2, -a) * Gamma(c)
        result /= Gamma(1 + 0.5*a - b) * Gamma(0.5 + 0.5*a)
        # If there has been an overflow, try again with logs
        if zisnan(result) or zisinf(result) or result == 0.0:
            result = exp(
                LOG_PI_2 + lgam(c) - a * M_LN2 +
                -lgam(1 + 0.5*a - b) - lgam(0.5 + 0.5*a)
            )
            result *= gammasgn(c) * gammasgn(1 + 0.5*a - b)
            result *= gammasgn(0.5 + 0.5*a)
        return result
    # Reduces to a polynomial when a or b is a negative integer.
    if a == floor(a) and a < 0:
        return hyp2f1_series(a, b, c, z, EPS, <int> fabs(a) + 1)
    if b == floor(b) and b < 0:
        return hyp2f1_series(a, b, c, z, EPS, <int> fabs(b) + 1)
    # If one of c - a or c - b is a negative integer, reduces to evaluating
    # a polynomial through an Euler hypergeometric transformation
    # (DLMF 15.8.1).
    # hyp2f1(a, b, c, z) = (1 - z)**(c - a - b)*hyp2f1(c - a, c - b, c , z)
    condition1 = c - a == floor(c - a) and c - a < 0
    condition2 = c - b == floor(c - b) and c - b < 0
    if condition1 or condition2:
        if condition1:
            max_iter = <int> fabs(c - a) + 1
        if condition2:
            max_iter = <int> fabs(c - b) + 1
        result = zpow(1 - z, c - a - b)
        result *= hyp2f1_series(c - a, c - b, c, z, EPS, max_iter)
        return result
    # |z| < 0, real(z)  >= 0. Use defining Taylor series.
    # --------------------------------------------------------------------------
    # Although the series is convergent in the interior of the unit disk,
    # other methods offer better convergence for |z| >= 0.9. The z/(z-1)
    # transform offers more rapid convergence for real(z) < 0.
    if modulus_z < 0.9 and z.real >= 0:
        # Apply Euler Hypergeometric Transformation (DLMF 15.8.1) to reduce
        # magnitude of a and b when this is possible. Calculation of hyp2f1
        # using series methods becomes inaccurate for larger |a| and |b| due
        # to slow convergence. In the original Fortran implementation from
        # specfun this step was incorrectly applied if c - a  < a and
        # c - b < b without taking  absolute values, hurting precision in some
        # cases.
        if fabs(c - a) < fabs(a) and fabs(c - b) < fabs(b):
            result = zpow(1 - z, c - a - b)
            # Maximum number of terms 1500 comes from Fortran original.
            result *= hyp2f1_series(c - a, c - b, c, z, EPS, 1500)
            return result
        # Maximum number of terms 1500 comes from Fortran original.
        return hyp2f1_series(a, b, c, z, EPS, 1500)
    # Fall through to original Fortran implementation.
    # -------------------------------------------------------------------------
    return double_complex_from_npy_cdouble(
        chyp2f1_wrap(a, b, c, npy_cdouble_from_double_complex(z))
    )


@cython.cdivision(True)
cdef inline double complex hyp2f1_series(
        double a,
        double b,
        double c,
        double complex z,
        double rtol,
        int max_iter
) nogil:
    cdef:
        int k
        double complex term = 1 + 0j
        double complex result = 1 + 0j
        double complex previous = 0 + 0j
    for k in range(max_iter):
        term *= z * (a + k) * (b + k) / ((c + k) * (k + 1))
        previous = result
        result += term
        if zabs(result - previous) <= rtol * zabs(result):
            break
    else:
        sf_error.error("hyp2f1", sf_error.NO_RESULT, NULL)
        result = zpack(NPY_NAN, NPY_NAN)
    return result
