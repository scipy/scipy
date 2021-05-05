# -*-cython-*-
#
# Implementation of the inverse of the logarithm of the CDF of the standard
# normal distribution.
#
# Copyright: Albert Steppi
#
# Distributed under the same license as SciPy
#
#
# Implementation Overview
#
# The inverse of the CDF of the standard normal distribution is available
# in scipy through the Cephes Math Library where it is called ndtri.
# We call our implementation of the inverse of the log CDF ndtri_exp.
# For -2 <= y <= log(1 - exp(-2)),  ndtri_exp is computed as ndtri(exp(y)).
# as ndtri(exp(y)).

# For 0 < p < exp(-2), the cephes implementation of ndtri uses an approximation
# for ndtri(p) which is a function of z = sqrt(-2.0 * log(p)). Letting
# y = log(p), for y < -2, ndtri_exp uses this approximation in log(p) directly.
# This allows this implementation to achieve high precision for very small y,
# whereas ndtri(exp(y)) becomes inf when exp(y) underflows for y < ~ -745.1.
#
# When p > 1 - exp(-2), the Cephes implementation of ndtri uses the symmetry
# of the normal distribution and calculates ndtri(p) as -ndtri(1 - p) allowing
# for the use of the same approximation. When y > log(1 - exp(-2)) this
# implementation calculates ndtri_exp as -ndtri(-exp1m(y)).
#
#
# Accuracy
#
# Cephes provides the following relative error estimates for ndtri
#                      Relative error:
# arithmetic   domain        # trials      peak         rms
#    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
#    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
#
# When y < -2, ndtri_exp must have relative error at least as small as the
# Cephes implementation of ndtri for p < exp(-2). It relies on the same
# approximation but does not have to lose precision by passing from p to log(p)
# before applying the approximation.
#
# Relative error of ndtri for values of the argument p near 1 can be much higher
# than claimed by the above chart. For p near 1, symmetry is exploited to
# replace the calculation of ndtri(p) with -ndtri(1 - p). The inverse of the
# normal CDF increases so rapidly near the endpoints of [0, 1] that the loss
# of precision incurred by the subtraction 1 - p due to limitations in binary
# approximation can make a significant difference in the results. Using
# version 9.3.0 targeting x86_64-linux-gnu we've observed the following
#
#                                               Estimated Relative Error
#  ndtri(1e-8)      = -5.612001244174789        ''
# -ndtri(1 - 1e-8)  = -5.612001243305505        1.55e-10
#  ndtri(1e-16)     = -8.222082216130435        ''
# -ndtri(1 - 1e-16) = -8.209536151601387        0.0015
#
# If expm1 is correctly rounded for y in [log(1 - exp(-2), 0), ndtri_exp(y)
# for y > log(1 - exp(-2)) should have the same relative error to ndtri(p)
# for p > 1 - exp(-2), though as seen above, even then this error may be higher
# than desired. IEEE-754 provides no guarantee on the accuracy of expm1
# however, therefore accuracy of ndtri_exp in this range is platform
# dependent.
#
# The case
#
# -2 <= y <= log(1 - exp(-2)) ~ -0.1454
#
# corresponds to ~ 0.135 <= p <= ~ 0.865
#
# The derivative of ndtri is
# sqrt(2 * pi) * exp(ndtri(x)**2 / 2). It is ~4.597 at x ~ 0.135, decreases
# monotonically to \sqrt(2 * pi) ~ 2.507 at x = 0 and increases monotonically
# again to ~4.597 at x ~ 0.865.
#
# It can be checked that all higher derivatives follow a similar pattern,
# their absolute value takes on a maximum (for this interval) at x ~ 0.135,
# decrease to a minimum at x = 0 and increases to the same maximum at x ~ 0.865.
# Derivatives of all orders are positive at x=log(1 - exp(-2)). Thus the worst
# possible loss of precision in ndtri(exp(x)) in the interval
# [0, log(1 - exp(-2))] due to error in calculating exp(x) must occur near
# x=log(1 - exp(-2)). By symmetry, the worst possible loss of precision in
# [-2, log(1 - exp(-2)] must occur near the endpoints. We may observe
# empirically that error at the endpoints due to exp is not substantial.
# Assuming that exp is accurate within +-ULP (unit of least precision),
# we've observed a value of at most ~6.0474e-16 for
#
# abs(ndtri(x + epsilon) - ndtri(x))
#
# for x near exp(-2) or 1 - exp(-2) and epsilon equal to the unit of least
# precision of x.

# (IEEE-754 provides no guarantee on the accuracy of exp, but for most
# compilers on most architectures an assumption of +-ULP should be
# reasonable.)
#
# The error here is on the order of the error in the Cephes implementation of
# ndtri itself, leading to an error profile that is still favorable.


import cython
from libc.float cimport DBL_MAX
from libc.math cimport exp, expm1, log, log1p, sqrt, HUGE_VAL

cdef extern from "cephes/polevl.h":
    double polevl(double x, const double coef[], int N) nogil
    double p1evl(double x, const double coef[], int N) nogil

from ._cephes cimport ndtri

@cython.cdivision(True)
cdef inline double _ndtri_exp_small_y(double y) nogil:
    """Return inverse of log CDF of normal distribution for very small y

    For p sufficiently small, the inverse of the CDF of the normal
    distribution can be approximated to high precision as a rational function
    in sqrt(-2.0 * log(p)).
    """
    cdef:
        double P1[9]
        double Q1[8]
        double P2[9]
        double Q2[8]
        double x, x0, x1, z
    P1[0] = 4.05544892305962419923
    P1[1] = 3.15251094599893866154e1
    P1[2] = 5.71628192246421288162e1
    P1[3] = 4.40805073893200834700e1
    P1[4] = 1.46849561928858024014e1
    P1[5] = 2.18663306850790267539
    P1[6] = -1.40256079171354495875e-1
    P1[7] =-3.50424626827848203418e-2
    P1[8] = -8.57456785154685413611e-4
    Q1[0] = 1.57799883256466749731e1
    Q1[1] = 4.53907635128879210584e1
    Q1[2] = 4.13172038254672030440e1
    Q1[3] = 1.50425385692907503408e1
    Q1[4] = 2.50464946208309415979
    Q1[5] = -1.42182922854787788574e-1
    Q1[6] = -3.80806407691578277194e-2
    Q1[7] = -9.33259480895457427372e-4
    P2[0] = 3.23774891776946035970
    P2[1] = 6.91522889068984211695
    P2[2] = 3.93881025292474443415
    P2[3] = 1.33303460815807542389
    P2[4] = 2.01485389549179081538e-1
    P2[5] = 1.23716634817820021358e-2
    P2[6] = 3.01581553508235416007e-4
    P2[7] = 2.65806974686737550832e-6
    P2[8] = 6.23974539184983293730e-9
    Q2[0] = 6.02427039364742014255
    Q2[1] = 3.67983563856160859403
    Q2[2] = 1.37702099489081330271
    Q2[3] = 2.16236993594496635890e-1
    Q2[4] = 1.34204006088543189037e-2
    Q2[5] = 3.28014464682127739104e-4
    Q2[6] = 2.89247864745380683936e-6
    Q2[7] = 6.79019408009981274425e-9
    # sqrt(-2 * y) is faster and has more precision but overflows when
    # y < -DBL_MAX * 0.5
    if y >= -DBL_MAX * 0.5:
        x = sqrt(-2 * y)
    else:
        x = 1.4142135623730950 * sqrt(-y) # sqrt(-2 * y)
    x0 = x - log(x) / x
    z = 1 / x
    if x < 8.0:
        x1 = z * polevl(z, P1, 8) / p1evl(z, Q1, 8)
    else:
        x1 = z * polevl(z, P2, 8) / p1evl(z, Q2, 8)
    return x1 - x0


cdef inline double ndtri_exp(double y) nogil:
    """Return inverse of logarithm of Normal CDF evaluated at y."""
    if y < -DBL_MAX:
        return -HUGE_VAL
    elif y < - 2.0:
        return _ndtri_exp_small_y(y)
    elif y > -0.14541345786885906: # log1p(-exp(-2))
        return -ndtri(-expm1(y))
    else:
        return ndtri(exp(y))
