"""
Implementation of Gauss's hypergeometric function for complex values. This is
an effort to incrementally translate the Fortran implementation from specfun.f
into Cython so that it's easier to maintain and it's easier to correct defects
in the original implementation.  Computation of Gauss's hypergeometric function
involves handling a patchwork of special cases. The goal is to translate the
cases into Cython little by little, falling back to the Fortran implementation
for the cases that are yet to be handled. By default the Fortran original is
followed as closely as possible except for situations when an improvement is
obvious. Attempts are made to document the why behind certain decisions made by
the original implementation, with references to the NIST Digital Library of
Mathematical Functions [1] added where they are appropriate. The review paper
by Pearson et al [2] is an excellent resource for best practices for numerical
computation of hypergeometric functions, and the intent is to follow this paper
when making improvements to and correcting defects in the original
implementation.

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
[3] Raimundas Vidunas, "Degenerate Gauss Hypergeometric Functions",
    Kyushu Journal of Mathematics, 2007, Volume 61, Issue 1, Pages 109-135,
[4] López, J.L., Temme, N.M. New series expansions of the Gauss hypergeometric
    function. Adv Comput Math 39, 349-365 (2013).
    https://doi.org/10.1007/s10444-012-9283-y
"""

cimport cython
from numpy cimport npy_cdouble
from libc.stdint cimport uint64_t, UINT64_MAX
from libc.math cimport fabs, exp, M_LN2, M_PI, pow, trunc

from . cimport sf_error
from ._cephes cimport Gamma, gammasgn, lgam
from ._complexstuff cimport (
    double_complex_from_npy_cdouble,
    npy_cdouble_from_double_complex,
    zabs,
    zisinf,
    zisnan,
    zpack,
    zpow,
)


cdef extern from "numpy/npy_math.h":
    double NPY_NAN
    double NPY_INFINITY


cdef extern from 'specfun_wrappers.h':
    npy_cdouble chyp2f1_wrap(
        double,
        double,
        double,
        npy_cdouble,
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
        double modulus_z
        double max_degree
        double complex result
        bint a_neg_int, b_neg_int, c_non_pos_int
        bint c_minus_a_neg_int, c_minus_b_neg_int
    modulus_z = zabs(z)
    a_neg_int = a == trunc(a) and a < 0
    b_neg_int = b == trunc(b) and b < 0
    c_non_pos_int = c == trunc(c) and c <= 0
    c_minus_a_neg_int = c - a == trunc(c - a) and c - a < 0
    c_minus_b_neg_int = c - b == trunc(c - b) and c - b < 0
    # Special Cases
    # -------------------------------------------------------------------------
    # Takes constant value 1 when a = 0 or b = 0, even if c is a non-positive
    # integer. This follows mpmath.
    if a == 0 or b == 0:
        return 1.0 + 0.0j
    # Equals 1 when z is 0, unless c is 0.
    if modulus_z == 0:
        if c != 0:
            return 1.0 + 0.0j
        else:
            # Returning real part NAN and imaginary part 0 follows mpmath.
            return zpack(NPY_NAN, 0)
    # Diverges when c is a non-positive integer unless a is an integer with
    # c <= a <= 0 or b is an integer with c <= b <= 0, (or z equals 0 with
    # c != 0) Cases z = 0, a = 0, or b = 0 have already been handled. We follow
    # mpmath in handling the degenerate cases where any of a, b, c are
    # non-positive integers. See [3] for a treatment of degenerate cases.
    if c_non_pos_int and not (
            a_neg_int and c <= a < 0 or b_neg_int and c <= b < 0
    ):
        return NPY_INFINITY + 0.0j
    # Diverges as real(z) -> 1 when c < a + b.
    # Todo: Actually check for overflow instead of using a fixed tolerance for
    # all parameter combinations like in the Fortran original. This will have to
    # wait until all points in the neighborhood of z = 1 are handled in Cython.
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
    if zabs(z + 1) < EPS and fabs(1 + a - b - c) < EPS:
        # The computation below has been simplified through
        # Legendre duplication for the Gamma function (DLMF 5.5.5).
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
    # If a and b are both negative integers, we take care to terminate
    # the series at a or b of smaller magnitude. This is to ensure proper
    # handling of situations like a < c < b <= 0, a, b, c all non-positive
    # integers, where terminating at a would lead to a term of the form 0 / 0.
    if a_neg_int or b_neg_int:
        if a_neg_int and b_neg_int:
            max_degree = fabs(a) - 1 if a > b else fabs(b) - 1
        elif a_neg_int:
            max_degree = fabs(a) - 1
        else:
            max_degree = fabs(b) - 1
        if max_degree <= UINT64_MAX:
            # This cast is OK because we've ensured max_degree will fit into
            # an int.
            return hyp2f1_series(
                a, b, c, z, <uint64_t> max_degree, False, 0
            )
        else:
            sf_error.error("hyp2f1", sf_error.NO_RESULT, NULL)
            return zpack(NPY_NAN, NPY_NAN)
    # If one of c - a or c - b is a negative integer, reduces to evaluating
    # a polynomial through an Euler hypergeometric transformation.
    # (DLMF 15.8.1)
    if c_minus_a_neg_int or c_minus_b_neg_int:
        max_degree = (
            fabs(c - b) if c_minus_b_neg_int else fabs(c - a)
        )
        if max_degree <= UINT64_MAX:
            result = zpow(1 - z, c - a - b)
            result *= hyp2f1_series(
                c - a, c - b, c, z, <uint64_t> max_degree, False, 0
            )
            return result
        else:
            sf_error.error("hyp2f1", sf_error.NO_RESULT, NULL)
            return zpack(NPY_NAN, NPY_NAN)
    # |z| < 0, real(z) >= 0. Use defining Taylor series.
    # --------------------------------------------------------------------------
    if modulus_z < 0.9 and z.real >= 0:
        # Apply Euler Hypergeometric Transformation (DLMF 15.8.1) to reduce
        # size of a and b if possible. We follow the original Fortran
        # implementation although there is very likely a better heuristic to
        # determine when this transformation should be applied. As it stands,
        # it hurts precision in some cases.
        if c - a < a and c - b < b:
            result = zpow(1 - z, c - a - b)
            # Maximum number of terms 1500 comes from Fortran original.
            result *= hyp2f1_series(c - a, c - b, c, z, 1500, True, EPS)
            return result
        # Maximum number of terms 1500 comes from Fortran original.
        return hyp2f1_series(a, b, c, z, 1500, True, EPS)
    if 0.9 <= modulus_z < 1.1 and zabs(1 - z) >= 0.9 and z.real >= 0:
        # This condition for applying Euler Transformation (DLMF 15.8.1)
        # was determined empirically to work better for this case than that
        # used in the original Fortran implementation for |z| < 0.9,
        # real(z) >= 0.
        if (
                c - a <= a and c - b < b or
                c - a < a and c - b <= b
        ):
            result = zpow(1 - z, c - a - b)
            result *= hyp2f1_lopez_temme_series(c - a, c - b, c, z, 1500, EPS)
            return result
        return hyp2f1_lopez_temme_series(a, b, c, z, 1500, EPS)
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
        uint64_t max_degree,
        bint early_stop,
        double rtol,
) nogil:
    """Return Truncated Maclaurin series for hyp2f1.

    Series is convergent for |z| < 1 but is only practical for numerical
    computation when |z| < 0.9.

    Parameters
    ----------
    a : double
    b : double
    c : double
    z : double complex
    max_degree : uint64_t
        Maximum degree of terms before truncating.
    early_stop : bint
    rtol : double
        If early_stop is True, truncate when
            |current_sum - previous_sum| <= rtol * |current_sum|
        If max_degree is reached before the stopping condition is satisfied,
        return nan. If early_stop is False, sum all terms up until and
        including that with degree max_degree.

    Returns
    -------
    double complex
    """
    cdef:
        uint64_t k
        double complex term = 1 + 0j
        double complex previous = 0 + 0j
        double complex result = 1 + 0j
    for k in range(max_degree + 1):
        # Follows the Fortran original exactly. Using *= degrades precision.
        # Todo: Check generated assembly to see what's going on.
        term = term * (a + k) * (b + k) / ((k + 1) * (c + k)) * z
        previous = result
        result += term
        if early_stop and zabs(result - previous) < zabs(result) * rtol:
            break
    else:
        # Return nan if max_degree is exceeded without convergence and
        # early_stop has been set to True.
        if early_stop:
            sf_error.error("hyp2f1", sf_error.NO_RESULT, NULL)
            result = zpack(NPY_NAN, NPY_NAN)
    return result


cdef inline double complex hyp2f1_lopez_temme_series(
        double a,
        double b,
        double c,
        double complex z,
        int max_degree,
        double rtol,
) nogil:
    """Lopez-Temme Series for Gaussian hypergeometric function [4].

    Converges for all z with real(z) < 1, including in the regions surrounding
    the points exp(+- i*pi/3) that are not covered by any of the standard
    transformations.
    """
    cdef:
        int n
        double phi_previous, phi
        double complex prefactor, previous, Z, result
    prefactor = zpow(1 - 0.5 * z, -a)
    # phi(n, b, c) = hyp2f1(-n, b, c, 2). It is computed through a linear
    # recurrence of degree 2. phi and phi_previous below are the initial
    # conditions of this recurrence.
    phi, phi_previous = 1 - 2 * b / c, 1.0
    previous = 1 + 0j
    Z = a * z / (z - 2)
    result = previous + Z * phi
    for n in range(2, max_degree):
        phi, phi_previous = (
            ((n - 1) * phi_previous - (2 * b - c) * phi) / (c + n - 1), phi
        )
        Z = Z * (a + n - 1) * z / ((z - 2) * n)
        previous, result = result, result + Z * phi
        if zabs(result - previous) <= rtol * zabs(result):
            result = prefactor * result
            break
    else:
        sf_error.error("hyp2f1", sf_error.NO_RESULT, NULL)
        result = zpack(NPY_NAN, NPY_NAN)
    return result
