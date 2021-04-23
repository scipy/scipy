import cython
from libc.math cimport exp, expm1, log, log1p, sqrt

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
    x = sqrt(-2.0 * y)
    x0 = x - log(x) / x
    z = 1 / x
    if x < 8.0:
        x1 = z * polevl(z, P1, 8) / p1evl(z, Q1, 8)
    else:
        x1 = z * polevl(z, P2, 8) / p1evl(z, Q2, 8)
    return x1 - x0


cdef inline double ndtri_exp(double y) nogil:
    """Return inverse of logarithm of Normal CDF evaluated at y."""
    if y < - 2.0:
        return _ndtri_exp_small_y(y)
    elif y > log1p(-exp(-2)):
        return -ndtri(-expm1(y))
    else:
        return ndtri(exp(y))
