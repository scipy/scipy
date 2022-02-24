"""Implementation of Gauss's hypergeometric function for complex values.

This implementation is based on the Fortran implementation by Shanjie Zhang and
Jianming Jin included in specfun.f [1]_.  Computation of Gauss's hypergeometric
function involves handling a patchwork of special cases. Zhang and Jin's
implementation is being replaced case by case with this Cython implementation,
which falls back to their implementation for the cases that are yet to be
handled. By default the Zhang and Jin implementation has been followed as
closely as possible except for situations where an improvement was
obvious. We've attempted to document the reasons behind decisions made by Zhang
and Jin and to document the reasons for deviating from their implementation
when this has been done. References to the NIST Digital Library of Mathematical
Functions [2]_ have been added where they are appropriate. The review paper by
Pearson et al [3]_ is an excellent resource for best practices for numerical
computation of hypergeometric functions. We have followed this review paper
when making improvements to and correcting defects in Zhang and Jin's
implementation. When Pearson et al propose several competing alternatives for a
given case, we've used our best judgment to decide on the method to use.

Author: Albert Steppi

Distributed under the same license as Scipy.

References
----------
.. [1] S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996
.. [2] NIST Digital Library of Mathematical Functions. http://dlmf.nist.gov/,
       Release 1.1.1 of 2021-03-15. F. W. J. Olver, A. B. Olde Daalhuis,
       D. W. Lozier, B. I. Schneider, R. F. Boisvert, C. W. Clark, B. R. Miller,
       B. V. Saunders, H. S. Cohl, and M. A. McClain, eds.
.. [3] Pearson, J.W., Olver, S. & Porter, M.A.
       "Numerical methods for the computation of the confluent and Gauss
       hypergeometric functions."
       Numer Algor 74, 821-866 (2017). https://doi.org/10.1007/s11075-016-0173-0
.. [4] Raimundas Vidunas, "Degenerate Gauss Hypergeometric Functions",
       Kyushu Journal of Mathematics, 2007, Volume 61, Issue 1, Pages 109-135,
.. [5] López, J.L., Temme, N.M. New series expansions of the Gauss hypergeometric
       function. Adv Comput Math 39, 349-365 (2013).
       https://doi.org/10.1007/s10444-012-9283-y
"""

cimport cython
from numpy cimport npy_cdouble
from libc.stdint cimport uint64_t, UINT64_MAX
from libc.math cimport fabs, exp, M_LN2, M_PI, M_1_PI, pow, sin, trunc

from . cimport sf_error
from ._cephes cimport Gamma, gammasgn, lanczos_sum_expg_scaled, lgam
from ._complexstuff cimport (
    double_complex_from_npy_cdouble,
    npy_cdouble_from_double_complex,
    zabs,
    zisinf,
    zisnan,
    zpack,
    zpow,
)

cdef extern from "cephes/lanczos.h":
    double lanczos_g

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


# Small value from Zhang and Jin's Fortran implementation.
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
    # Note that in the following three conditions, we avoid the case where
    # c is a non-positive integer. This is because if c is a non-positive
    # integer and we have not returned already, then one of a or b must
    # be a negative integer and the desired result can be computed as a
    # polynomial with finite sum.
    if (
            fabs(1 - z.real) < EPS and z.imag == 0 and c - a - b < 0 and
            not c_non_pos_int
    ):
        return NPY_INFINITY + 0.0j
    # Gauss's Summation Theorem for z = 1; c - a - b > 0 (DLMF 15.4.20).
    if z == 1.0 and c - a - b > 0 and not c_non_pos_int:
        return gamma_ratio(c, c - a - b, c - a, c - b)
    # Kummer's Theorem for z = -1; c = 1 + a - b (DLMF 15.4.26).
    if zabs(z + 1) < EPS and fabs(1 + a - b - c) < EPS and not c_non_pos_int:
        return gamma_ratio(
            a - b + 1, 0.5*a + 1, a + 1, 0.5*a - b + 1
        )
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
        # size of a and b if possible. We follow Zhang and Jin's
        # implementation [1] although there is very likely a better heuristic
        # to determine when this transformation should be applied. As it
        # stands, it hurts precision in some cases.
        if c - a < a and c - b < b:
            result = zpow(1 - z, c - a - b)
            # Maximum number of terms 1500 comes from Zhang and Jin.
            result *= hyp2f1_series(c - a, c - b, c, z, 1500, True, EPS)
            return result
        # Maximum number of terms 1500 comes from Zhang and Jin.
        return hyp2f1_series(a, b, c, z, 1500, True, EPS)
    # Points near exp(iπ/3), exp(-iπ/3) not handled by any of the standard
    # transformations. Use series of López and Temme [5]. These regions
    # were not correctly handled by Zhang and Jin's implementation.
    # -------------------------------------------------------------------------
    if 0.9 <= modulus_z < 1.1 and zabs(1 - z) >= 0.9 and z.real >= 0:
        # This condition for applying Euler Transformation (DLMF 15.8.1)
        # was determined empirically to work better for this case than that
        # used in Zhang and Jin's implementation for |z| < 0.9,
        # real(z) >= 0.
        if (
                c - a <= a and c - b < b or
                c - a < a and c - b <= b
        ):
            result = zpow(1 - z, c - a - b)
            result *= hyp2f1_lopez_temme_series(c - a, c - b, c, z, 1500, EPS)
            return result
        return hyp2f1_lopez_temme_series(a, b, c, z, 1500, EPS)
    # z/(z - 1) transformation (DLMF 15.8.1). Avoids cancellation issues that
    # occur with Maclaurin series for real(z) < 0.
    # ------------------------------------------------------------------------
    if modulus_z < 1.1 and z.real < 0:
        if 0 < b < a < c:
            a, b = b, a
        return zpow(1 - z, -a) * hyp2f1_series(
            a, c - b, c, z/(z - 1), 500, True, EPS
        )
    # Fall through to Zhang and Jin's Fortran implementation.
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
        # Follows Zhang and Jin [1] exactly. Using *= degrades precision.
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


cdef inline double gamma_ratio(
        double u, double v, double w, double x
) nogil:
    """Computes a ratio of gamma functions.

    Computes gamma(u)*gamma(v)/(gamma(w)*gamma(x)).

    It is assumed that x = u + v - w, but it is left to the user to
    ensure this. Calculating u + v - w within the function itself can
    lead to unnecessary floating point error.

    Checks for overflow or underflow and trys again using the lanczos
    approximation if needed.
    """
    cdef double result
    result = Gamma(u) * Gamma(v)
    result /= Gamma(w) * Gamma(x)
    # Try again with lanczos approximation if there has been
    # underflow or overflow.
    if zisnan(result) or zisinf(result) or result == 0.0:
        result = gamma_ratio_lanczos(u, v, w, x)
    return result
    

@cython.boundscheck(False)
@cython.cdivision(True)
cdef inline double gamma_ratio_lanczos(
        double u, double v, double w, double x
) nogil:
    """Compute ratio of gamma functions using lanczos approximation.

    Computes gamma(u)*gamma(v)/(gamma(w)*gamma(x))

    It is assumed that x = u + v - w, but it is left to the user to
    ensure this. Calculating u + v - w within the function itself can
    lead to unnecessary floating point error.

    The lanczos approximation takes the form

    gamma(x) = fac * lanczos_sum_expg_scaled(x)
    
    where fac = ((x + lanczos_g - 0.5)/e)**(x - 0.5).

    The terms can be combined analytically to avoid underflow and overflow.
    The formula above is only valid for x >= 0.5, but can be extended to
    x < 0.5 with the reflection principle. This implementation priortizes
    simplifying the logic of when to apply the reflection principle over
    using the most efficient analytical combination of terms for each
    case.
    """
    cdef:
        double g
        double lanczos_part,
        double factor_part
        double factors[4]  # Lanczos factors for gamma u, v, w, x respectively
        double min_factor
        double max_factor
        double factors_midpoint
        double current_distance_to_midpoint
        double min_distance_to_midpoint
        int i
        int absorbed_index
    # The below implementation may incorrectly return finite results
    # at poles of the gamma function. Handle these cases explicitly.
    if u == trunc(u) and u <= 0 or v == trunc(v) and v <= 0:
        # Return nan if numerator has pole. Diverges to +- infinity
        # depending on direction so value is undefined.
        return NPY_NAN
    if w == trunc(w) and w <= 0 or x == trunc(x) and x <= 0:
        # Return 0 if denominator has pole but not numerator.
        return 0
    g = lanczos_g
    lanczos_part = 1
    factor_part = 1
    # Due to the argument of the second gamma function in the denominator
    # being u + v - w, the factors of the lanczos approximation cancel
    # nicely. The impact of applying the reflection principle when an
    # argument is less than 0.5 is relatively straightforward, but the
    # logic of the code is somewhat tricky. The best way to understand is
    # to study the code from back to front, starting with the result and
    # working backwards through how it was calculated.
    if u >= 0.5:
        lanczos_part *= lanczos_sum_expg_scaled(u)
        factors[0] = (u + g - 0.5)
    else:
        lanczos_part /= lanczos_sum_expg_scaled(1 - u)*sin(M_PI*u)*M_1_PI
        factors[0] = (1 - u + g - 0.5)
    if v >= 0.5:
        lanczos_part *= lanczos_sum_expg_scaled(v)
        factors[1] = (v + g - 0.5)
    else:
        lanczos_part /= lanczos_sum_expg_scaled(1 - v)*sin(M_PI*v)*M_1_PI
        factors[1] = (1 - v + g - 0.5)
    if w >= 0.5:
        lanczos_part /= lanczos_sum_expg_scaled(w)
        factors[2] = (w + g - 0.5)
    else:
        lanczos_part *= lanczos_sum_expg_scaled(1 - w)*sin(M_PI*w)*M_1_PI
        factors[2] = (1 - w + g - 0.5)
    if x >= 0.5:
        lanczos_part /= lanczos_sum_expg_scaled(x)
        factors[3] = (x + g - 0.5)
    else:
        lanczos_part *= lanczos_sum_expg_scaled(1 - x)
        lanczos_part *= sin(M_PI*x)*M_1_PI
        factors[3] = (1 - x + g - 0.5)
    # Decide on how to combine terms by finding factor closest to
    # midpoint between the largest and smallest factors.
    min_factor = NPY_INFINITY
    max_factor = -NPY_INFINITY
    for i in range(4):
        if factors[i] < min_factor:
            min_factor = factors[i]
        if factors[i] > max_factor:
            max_factor = factors[i]
    # Midpoint calculation uses extra floating point op to guard against
    # overflow
    factors_midpoint = min_factor / 2 + max_factor / 2
    min_distance_to_midpoint = NPY_INFINITY
    # Identify factor closest to the midpoint
    for i in range(4):
        current_distance_to_midpoint = fabs(factors[i] - factors_midpoint)
        if (
                current_distance_to_midpoint < min_distance_to_midpoint or
                current_distance_to_midpoint == min_distance_to_midpoint and
                factors[i] > factors[absorbed_index]
        ):
            min_distance_to_midpoint = current_distance_to_midpoint
            absorbed_index = i
    # Absorb u factor into others.
    if absorbed_index == 0:
        factor_part *= pow(factors[1] / factors[0], v - 0.5)
        factor_part *= pow(factors[0] / factors[2], w - 0.5)
        factor_part *= pow(factors[0] / factors[3], x - 0.5)
    # Absorb v factor into others.
    elif absorbed_index == 1:
        factor_part *= pow(factors[0] / factors[1], u - 0.5)
        factor_part *= pow(factors[1] / factors[2], w - 0.5)
        factor_part *= pow(factors[1] / factors[3], x - 0.5)
    # Absorb w factor into others.
    elif absorbed_index == 2:
        factor_part *= pow(factors[0] / factors[2], u - 0.5)
        factor_part *= pow(factors[1] / factors[2], v - 0.5)
        factor_part *= pow(factors[2] / factors[3], x - 0.5)
    # Absorb x factor into others.
    elif absorbed_index == 3:
        factor_part *= pow(factors[0]/factors[3], u - 0.5)
        factor_part *= pow(factors[1]/factors[3], v - 0.5)
        factor_part *= pow(factors[3]/factors[2], w - 0.5)
    return factor_part * lanczos_part
