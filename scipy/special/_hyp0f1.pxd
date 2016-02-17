from libc.math cimport pow, sqrt, floor, log, log1p, exp, M_PI, fabs
from numpy.math cimport NAN, isinf
cimport numpy as np

from _xlogy cimport xlogy
from _complexstuff cimport zsqrt, zpow, zabs

cdef extern from "float.h":
    double DBL_MAX, DBL_MIN

cdef extern from "cephes.h":
    double iv(double v, double x) nogil
    double jv(double n, double x) nogil
    double Gamma(double x) nogil
    double lgam(double x) nogil

cdef extern from "c_misc/misc.h":
    double gammasgn(double x) nogil

cdef extern from "amos_wrappers.h":
    np.npy_cdouble cbesi_wrap(double v, np.npy_cdouble z) nogil
    np.npy_cdouble cbesj_wrap(double v, np.npy_cdouble z) nogil
    double sin_pi(double x) nogil

#
# Real-valued kernel
#
cdef inline double _hyp0f1_real(double v, double z) nogil:
    cdef double arg, v1, arg_exp, bess_val

    # handle poles, zeros 
    if v <= 0.0 and v == floor(v):
        return NAN
    if z == 0.0 and v != 0.0:
        return 1.0

    # both v and z small: truncate the Taylor series at O(z**2)
    if fabs(z) < 1e-6*(1.0 + fabs(v)):
        return 1.0 + z/v + z*z/(2.0*v*(v+1.0))

    if z > 0:
        arg = sqrt(z)
        arg_exp = xlogy(1.0-v, arg) + lgam(v)
        bess_val = iv(v-1, 2.0*arg)
        
        if (arg_exp > log(DBL_MAX) or bess_val == 0 or   # overflow
            arg_exp < log(DBL_MIN) or isinf(bess_val)):  # underflow
            return _hyp0f1_asy(v, z)
        else:
            return exp(arg_exp) * gammasgn(v) * bess_val
    else:
        arg = sqrt(-z)
        return pow(arg, 1.0 - v) * Gamma(v) * jv(v - 1, 2*arg)


cdef inline double _hyp0f1_asy(double v, double z) nogil:
    r"""Asymptotic expansion for I_{v-1}(2*sqrt(z)) * Gamma(v) 
    for real $z > 0$ and $v\to +\infty$.

    Based off DLMF 10.41
    """
    cdef:
        double arg = sqrt(z)
        double v1 = fabs(v - 1)
        double x = 2.0 * arg / v1
        double p1 = sqrt(1.0 + x*x)
        double eta = p1 + log(x) - log1p(p1)
        double arg_exp_i, arg_exp_k
        double pp, p2, p4, p6, u1, u2, u3, u_corr_i, u_corr_k
        double result, gs

    arg_exp_i = -0.5*log(p1)
    arg_exp_i -= 0.5*log(2.0*M_PI*v1)
    arg_exp_i += lgam(v)
    gs = gammasgn(v)

    arg_exp_k = arg_exp_i
    arg_exp_i += v1 * eta
    arg_exp_k -= v1 * eta

    # large-v asymptotic correction, DLMF 10.41.10
    pp = 1.0/p1
    p2 = pp*pp
    p4 = p2*p2
    p6 = p4*p2
    u1 = (3.0 - 5.0*p2) * pp / 24.0
    u2 = (81.0 - 462.0*p2 + 385.0*p4) * p2 / 1152.0
    u3 = (30375.0 - 369603.0*p2 + 765765.0*p4 - 425425.0*p6) * pp * p2 / 414720.0
    u_corr_i = 1.0 + u1/v1 + u2/(v1*v1) + u3/(v1*v1*v1)

    result = exp(arg_exp_i - xlogy(v1, arg)) * gs * u_corr_i
    if v - 1 < 0:
        # DLMF 10.27.2: I_{-v} = I_{v} + (2/pi) sin(pi*v) K_v
        u_corr_k = 1.0 - u1/v1 + u2/(v1*v1) - u3/(v1*v1*v1)
        result += exp(arg_exp_k + xlogy(v1, arg)) * gs * 2.0 * sin_pi(v1) * u_corr_k

    return result


#
# Complex valued kernel
#
cdef inline double complex _hyp0f1_cmplx(double v, double complex z) nogil:
    cdef:
        np.npy_cdouble zz = (<np.npy_cdouble*>&z)[0]
        np.npy_cdouble r
        double complex arg, s

    # handle poles, zeros 
    if v <= 0.0 and v == floor(v):
        return NAN
    if z.real == 0.0 and z.imag == 0.0 and v != 0.0:
        return 1.0

    # both v and z small: truncate the Taylor series at O(z**2)
    if zabs(z) < 1e-6*(1.0 + fabs(v)):
        return 1.0 + z/v + z*z/(2.0*v*(v+1.0))

    if zz.real > 0:
        arg = zsqrt(z)
        s = 2.0 * arg
        r = cbesi_wrap(v-1.0, (<np.npy_cdouble*>&s)[0] )
    else:
        arg = zsqrt(-z)
        s = 2.0 * arg
        r = cbesj_wrap(v-1.0, (<np.npy_cdouble*>&s)[0] )

    return (<double complex*>&r)[0] * Gamma(v) * zpow(arg, 1.0 - v)

