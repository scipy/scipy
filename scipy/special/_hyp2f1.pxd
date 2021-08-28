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
"""

cimport cython
from numpy cimport npy_cdouble
from libc.math cimport fabs, trunc

from . cimport sf_error
from ._complexstuff cimport (
    double_complex_from_npy_cdouble,
    npy_cdouble_from_double_complex,
    zabs,
    zpack,
)


cdef extern from "limits.h":
    cdef int INT_MAX


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


@cython.cdivision(True)
cdef inline double complex hyp2f1_complex(
        double a, double b, double c, double complex z
) nogil:
    cdef:
        int num_terms
        double modulus_z
        bint a_neg_int, b_neg_int, c_non_pos_int
    modulus_z = zabs(z)
    a_neg_int = a == trunc(a) and a < 0
    b_neg_int = b == trunc(b) and b < 0
    c_non_pos_int = c == trunc(c) and c <= 0
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
    # Fallback to Fortran original. To be translated to Cython later.
    if z == 1.0 and c - a - b > 0:
        return double_complex_from_npy_cdouble(
            chyp2f1_wrap(a, b, c, npy_cdouble_from_double_complex(z))
        )
    # Kummer's Theorem for z = -1; c = 1 + a - b (DLMF 15.4.26).
    # Fall back to Fortran original. To be translated to Cython later.
    if zabs(z + 1) < EPS and fabs(1 + a - b - c) < EPS:
        return double_complex_from_npy_cdouble(
            chyp2f1_wrap(a, b, c, npy_cdouble_from_double_complex(z))
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
        if max_degree <= INT_MAX:
            # This cast is OK because we've ensured max_degree will fit into
            # an int.
            return hyp2f1_series(
                a, b, c, z, <int> max_degree, False, 0
            )
        else:
            sf_error.error("hyp2f1", sf_error.NO_RESULT, NULL)
            return zpack(NPY_NAN, NPY_NAN)
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
        int max_degree,
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
    max_degree : int
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
        int k
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
