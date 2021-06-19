"""
Implementation of Gauss's hypergeometric function for complex values. This is
an effort to incrementally translate the Fortran implementation from specfun.f
into Cython so that it's easier to maintain and it's easier to correct defects
in the original implementation.

Computation of Gauss's hypergeometric function involves handling a patchwork of
special cases. The goal is to translate the cases into Cython little by little,
falling back to the Fortran implementation for the cases that are yet to be
handled. By default the Fortran original is followed as closely as possible
except for situations when an improvement is obvious. Attempts are made to
document the why behind certain decisions made by the original implementation,
with references to the NIST Digital Library of Mathematical Functions [1] added
where they are appropriate. The review paper by Pearson et al [2] is an
excellent resource for best practices for numerical computation of
hypergeometric functions, and the intent is to follow this paper when making
improvements to and correcting defects in the original implementation.

Author: Albert Steppi

Distributed under the same license as Scipy.

References
----------
[1] NIST Digital Library of Mathematical Functions. http://dlmf.nist.gov/,
    Release 1.1.1 of 2021-03-15. F. W. J. Olver, A. B. Olde Daalhuis,
    D. W. Lozier, B. I. Schneider, R. F. Boisvert, C. W. Clark, B. R. Miller,
    B. V. Saunders, H. S. Cohl, and M. A. McClain, eds.

[2] Pearson, J.W., Olver, S. & Porter, M.A.
    "Numerical methods for the computation of the confluent and Gauss
    hypergeometric functions."
    Numer Algor 74, 821-866 (2017). https://doi.org/10.1007/s11075-016-0173-0
"""


cimport cython
from numpy cimport npy_cdouble
from libc.math cimport fabs, exp, M_LN2, M_PI, pow, trunc

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

# Small value from the Fortran original.
DEF EPS = 1e-15

DEF SQRT_PI = 1.7724538509055159  # sqrt(M_PI)
DEF LOG_PI_2 = 0.5723649429247001  # log(M_PI) / 2


@cython.cdivision(True)
cdef inline double complex hyp2f1_complex(
        double a, double b, double c, double complex z
) nogil:
    cdef:
        int num_terms
        double modulus_z
        double complex result
        bint condition1, condition2
    modulus_z = zabs(z)
    # Special Cases
    # -------------------------------------------------------------------------
    # Equals 1 at z = 0. Takes constant value 1 when a = 0 or b = 0.
    if modulus_z == 0 or a == 0 or b == 0:
        return 1.0 + 0.0j
    # Diverges when c is a negative integer unless c <= a <= 0 or c <= b <= 0
    # or z = 0. Cases z = 0, a = 0, or b = 0 have already been handled. The
    # original Fortran implementation did not handle the case where c is a
    # negative integer correctly, returning infinity in all situations.
    if c == trunc(c) and c < 0 and not ((a == trunc(a) and c <= a < 0)
                                        or (b == trunc(b) and c <= b < 0)):
        return NPY_INFINITY + 0.0j
    # Diverges as real(z) -> 1 when c < a + b.
    if fabs(1 - z.real) < EPS and z.imag == 0 and c - a - b < 0:
        return NPY_INFINITY + 0.0j
    # Gauss's Summation Theorem for z = 1; c - a - b > 0 (DLMF 15.4.20).
    if z == 1.0 and c - a - b > 0:
        result = Gamma(c) * Gamma(c - a - b)
        result /= Gamma(c - a) * Gamma(c - b)
        # Try again with logs if there has been an overflow.
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
        # Try again with logs if there has been an overflow.
        if zisnan(result) or zisinf(result) or result == 0.0:
            result = exp(
                LOG_PI_2 + lgam(c) - a * M_LN2 +
                -lgam(1 + 0.5*a - b) - lgam(0.5 + 0.5*a)
            )
            result *= gammasgn(c) * gammasgn(1 + 0.5*a - b)
            result *= gammasgn(0.5 + 0.5*a)
        return result
    # Reduces to a polynomial when a or b is a negative integer.
    condition1 = a == trunc(a) and a < 0
    condition2 = b == trunc(b) and b < 0
    # This has become very tricky. If c is a negative integer but c == a,
    # then the series converges but if we naively calculate the series up
    # to the term of degree |a|, we get a 0 / 0 for the last term. We need
    # to stop at the term of degree |a| - 1. A similar consideration applies
    # if c == b. If a, b, and c are all negative integers with one of a or
    # b of smaller magnitude than c, we need to pick a or b of smallest
    # magnitude to determine the stopping point.
    if condition1 or condition2:
        if (condition1 and condition2):
            if a > b:
                num_terms = <int> fabs(a) - 1
            else:
                num_terms = <int> fabs(b) - 1
        elif condition1:
            num_terms = <int> fabs(a) - 1
        else:
            num_terms = <int> fabs(b) - 1
        return hyp2f1_series_fixed(a, b, c, z, num_terms)
    # If one of c - a or c - b is a negative integer, reduces to evaluating
    # a polynomial through an Euler hypergeometric transformation
    # (DLMF 15.8.1).
    # hyp2f1(a, b, c, z) = (1 - z)**(c - a - b)*hyp2f1(c - a, c - b, c , z)
    condition1 = c - a == trunc(c - a) and c - a < 0
    condition2 = c - b == trunc(c - b) and c - b < 0
    if condition1 or condition2:
        if condition1:
            num_terms = <int> fabs(c - a)
        if condition2:
            num_terms = <int> fabs(c - b)
        result = zpow(1 - z, c - a - b)
        result *= hyp2f1_series_fixed(c - a, c - b, c, z, num_terms)
        return result
    # |z| < 0, real(z) >= 0. Use defining Taylor series.
    # --------------------------------------------------------------------------
    # Although the series is convergent in the interior of the unit disk,
    # other methods offer better convergence for |z| >= 0.9. The z/(z-1)
    # transform offers better convergence for real(z) < 0.
    if modulus_z < 0.9 and z.real >= 0:
        # Apply Euler Hypergeometric Transformation (DLMF 15.8.1) to reduce
        # size of a and b if possible. Convergence slows as a and b grow in
        # the positive direction. For negative a or b, things are trickier and
        # its possible that this transformation can hurt precision by making
        # a or b take negative values that are large in magnitude. We follow
        # the Fortran original rather than trying to work out a better
        # heuristic for when this transformation improves precision.
        if c - a < a and c - b < b:
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
    """Maclaurin Series for hyp2f1.

    Stops computing when the modulus of the difference between the current
    value of the series and the previous value is smaller than tolerance rtol
    times the modulus of the current value. This stopping criterion is copied
    from the Fortran original. Returns nan if max_iter is reached before
    converging.
    """
    cdef:
        int k
        double complex term = 1 + 0j
        double complex result = 1 + 0j
        double complex previous
    for k in range(max_iter + 1):
        previous = result
        # Follows the Fortran original exactly. Do not rewrite with
        # term *=, this degrades performance with GCC 10.2.0 running on
        # x86_64-linux-gnu.
        term = term * (a + k) * (b + k) / ((k + 1) * (c + k)) * z
        result += term
        if zabs(result - previous) <= zabs(result) * rtol:
            break
    else:
        sf_error.error("hyp2f1", sf_error.NO_RESULT, NULL)
        result = zpack(NPY_NAN, NPY_NAN)
    return result


@cython.cdivision(True)
cdef inline double complex hyp2f1_series_fixed(
        double a,
        double b,
        double c,
        double complex z,
        int num_terms
) nogil:
    """Maclaurin series for hyp2f1 truncated after a fixed number of terms.

    Used when any of a, b, c - a, or c - a are negative integers and the series
    reduces to a polynomial. This allows this implementation to follow the
    Fortran original more closely. The original does not use early stopping
    when a tolerance is reached for these cases.
    """
    cdef:
        int k
        double complex term = 1 + 0j
        double complex result = 1 + 0j
    for k in range(num_terms + 1):
        # Follows the Fortran original exactly. Do not rewrite with
        # term *=, this degrades performance with GCC 10.2.0 running on
        # x86_64-linux-gnu.
        term = term * (a + k) * (b + k) / ((k + 1) * (c + k)) * z
        result += term
    return result
