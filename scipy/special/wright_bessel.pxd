# -*-cython-*-
#
# Implementation of Wright's generalized Bessel function Phi [1].
#
# Copyright: Christian Lorentzen
#
# Distributed under the same license as SciPy
#
#
# Implementation Overview:
#
# So far, only non-negative values of rho=a, beta=b and z=x are implemented.
# There are 4 different approaches depending on the ranges of the arguments.
# 1. Tylor series expansion in x=0 [1]. Suited for x <= 1. Involves gamma
#    funtion in each term.
# 2. Taylor series in a=0. Suited for tiny a and not too large x.
# 3. Asymptotic expansion for large x [2, 3]. Suited for large x while still
#    small a and b.
# 4. Integral representation [4]. Suited for all arguments.
#
# References:
# [1] https://dlmf.nist.gov/10.46.E1
# [2] E. M. Wright (1935), The asymptotic expansion of the generalized Bessel
#     function. Proc. London Math. Soc. (2) 38, pp. 257–270.
#     https://doi.org/10.1112/plms/s2-38.1.257
# [3] R. B. Paris (2017), The asymptotics of the generalised Bessel function,
#     Mathematica Aeterna, Vol. 7, 2017, no. 4, 381 - 406,
#     https://arxiv.org/abs/1711.03006
# [4] Y. F. Luchko (2008), Algorithms for Evaluation of the Wright Function for
#     the Real Arguments’ Values, Fractional Calculus and Applied Analysis 11(1)
#     http://sci-gems.math.bas.bg/jspui/bitstream/10525/1298/1/fcaa-vol11-num1-2008-57p-75p.pdf

import cython
from libc.math cimport cos, exp, floor, log10, pow, sin, sqrt, M_PI

cdef extern from "_c99compat.h":
    int sc_isnan(double x) nogil
    int sc_isinf(double x) nogil

from ._cephes cimport Gamma, rgamma, zeta, sinpi, cospi
from ._digamma cimport digamma
from ._complexstuff cimport inf, nan
from . cimport sf_error

# rgamma_zero: smallest value x for which rgamma(x) == 0 as x gets large
DEF rgamma_zero = 178.47241115886637
# exp_inf: smallest value x for which exp(x) == inf
DEF exp_inf = 709.78271289338403
# exp_zero: largest value x for which exp(x) == 0
DEF exp_zero = -745.13321910194117


# 1. Tylor series expansion in x=0
#    Phi(a, b, x) = sum_k x^k / k! / Gamma(a*k+b)
#
# Note that every term, and therefore also Phi(a, b, x), is monotone decreasing
# with increasing a or b.
@cython.cdivision(True)
cdef inline double _wb_series(double a, double b, double x,
    unsigned int nstart, unsigned int nstop) nogil:
    cdef:
        unsigned int k, k_max
        double xk_k, res

    xk_k = x**nstart * Gamma(nstart + 1)  # x^k/k!
    res = xk_k * rgamma(nstart*a + b)
    # term k=nstart+1 , +2 +3, ...
    if nstop > nstart:
        # series expansion until term k such that a*k+b <= rgamma_zero
        k_max = int(floor((rgamma_zero - b) / a))
        if nstop > k_max:
            nstop = k_max
        for k in range(nstart + 1, nstop):
            xk_k *= x/k
            res += xk_k * rgamma(a*k + b)

    return res


# 2. Tylor series expansion in a=0 up to order 2
#
# import sympy
# from sympy import S, Sum, factorial, gamma, symbols
# a, b, z, k = symbols("a b z k")
# expression = Sum(z**k/factorial(k)/gamma(a*k+b), (k, 0, S.Infinity))
# # nth term
# n = 1
# expression.diff(a, n).subs(a, 0).simplify().doit()
#
#    Phi(a, b, x) = exp(z)/Gamma(b)
#                   * (1 - a*z*Psi(b) + a^2/2*z(1+z)(Psi(b)^2 - Psi'(b))
#                      + O(a^3))
#    where Psi is the digamma function.
@cython.cdivision(True)
cdef inline double _wb_small_a(double a, double b, double x, int order) nogil:
    cdef double dg, dg1, dg2, dg3, res
    dg = digamma(b)
    # dg1 = polygamma(1, b)
    dg1 = zeta(2, x)
    # res = 1 - a*x*dg + a**2/2*x*(1+x)*(dg**2 - dg1)
    res = 1 + a*x*(-dg + 0.5*a*(1 + x)*(dg**2 - dg1))
    if order >= 3:
        # dg2 = polygamma(2, b)
        dg2 = -2. * zeta(3, x)
        res += -a**3/6. * x*(x**2 + 3*x + 1)*(dg**3 - 3*dg*dg1 + dg2)
    if order >= 4:
        # dg3 = polygamma(3, b)
        dg3 =  6 * zeta(4, x)
        res += a**4/24. * x*(x**3 + 6*x**2 + 7*x + 1) \
               *(dg**4 - 6*dg**2*dg1 + 4*dg*dg2 + 3*dg1**2 - dg3)
    res *= exp(x) * rgamma(b)
    return res


# 3. Asymptotic expansion for large x
#    Phi(a, b, x) ~ Z^(1/2-b) * exp((1+a)/a * Z) * sum_k (-1)^k * a_k / Z^k
#    with Z = (a*x)**(1/(1+a)).
#
# Wright (1935) lists the coefficients a_0 and a_1. With slightly different
# notation, Paris (2017) lists coefficients c_k up to order k=3.
# Paris (2017) uses ZP = (1+a)/a * Z  (ZP = Z of Paris) and
# a_k = a_0 * (-a/(1+a))^k c_k
#
# The coefficients a_k can be generated with the following code:
#
# import sympy
# from sympy import factorial, gamma, gammasimp, pi, symbols, Function, Rational
# class g(Function):
#     nargs = 3
#
#     @classmethod
#     def eval(cls, n, rho, v):
#         if not n >= 0:
#             raise ValueError("must have n >= 0")
#
#         if n == 0:
#             return 1
#         else:
#             return g(n-1, rho,v) + gammasimp(gamma(rho+2+n)/gamma(rho+2)) \
#                                  / gammasimp(gamma(3+n)/gamma(3))*v**n
#
# class a(Function):
#     nargs = 3
#
#     @classmethod
#     def eval(cls, m, rho, beta):
#         if not m >= 0:
#             raise ValueError("must have m >= 0")
#
#         v = symbols("v")
#         expression = (1-v)**(-beta) * g(2*m, rho, v)**(-m-Rational(1, 2))
#         res = expression.diff(v, 2*m).subs(v, 0) / factorial(2*m)
#         res = res * (gamma(m + Rational(1, 2)) / (2*pi) \
#                      * (2/(rho+1))**(m + Rational(1, 2)))
#         return res
#
# m, rho, beta = symbols("m rho beta")
# a0 = a(0, rho, beta)
# a1 = a(1, rho, beta)
#
@cython.cdivision(True)
cdef inline double _wb_asymptotic(double a, double b, double x) nogil:
    cdef:
        double A[9]  # powers of a
        double B[11]  # powers of b
        double Ap1[6]  # powers of (1+a)
        double C[6]  # coefficients of asymptotic series a_k
        double Z, Zp, res
        int k

    A[0] = 1.
    B[0] = 1.
    Ap1[0] = 1.
    for k in range(1, 9):
        A[k] = A[k-1] * a
    for k in range(1, 11):
        B[k] = B[k-1] * b
    for k in range(1, 6):
        Ap1[k] = Ap1[k-1] * (1 + a)

    C[0] = 1./sqrt(2. * M_PI * Ap1[1])

    C[1] = C[0] / (24 * Ap1[1])
    C[1] *= (2*a + 1)*(2 + a) - 12*b*(1 + a - b)

    C[2] = C[0] / (1152 * Ap1[2])
    C[2] *= ((2 + a)*(1 + 2*a)*(2 - 19*a + 2*A[2]) \
            + 24*b*Ap1[1]*(2 + 7*a - 6*A[2]) \
            - 24*B[2]*(4 - 5*a - 20*A[2]) - 96*B[3]*(1 + 5*a) + 144*B[4])

    C[3] = -C[0] / (414720*Ap1[3])
    C[3] *= (2 + a)*(1 + 2*a)*(556 + 1628*a - 9093*A[2] + 1628*A[3] \
                               + 556*A[4]) \
            - 180*b*Ap1[1]*(12 - 172*a - 417*A[2] + 516*A[3] - 20*A[4]) \
            - 180*B[2]*(76 + 392*a - 567*A[2] - 1288*A[3] + 364*A[4]) \
            + 1440*B[3]*(8 - 63*a - 147*A[2] + 112*A[3]) \
            + 10800*B[4]*(2 + 7*a - 14*A[2]) \
            - 8640*B[5]*(1 - 7*a) - 8640*B[6]

    C[4] = C[0] / (39813120*Ap1[4])
    C[4] *= 103680*B[8] - 414720*B[7]*(3*a - 1) + 725760*B[6]*a*(8*a - 7) \
            - 48384*B[5]*(274*A[3] - 489*A[2] + 39*a + 26) \
            + 30240*B[4]*(500*A[4] - 1740*A[3] + 495*A[2] + 340*a - 12) \
            - 2880*B[3]*(2588*A[5] - 19780*A[4] + 14453*A[3] + 9697*A[2] \
                         - 1892*a - 404) \
            + 48*B[2]*(11488*A[6] - 547836*A[5] + 1007484*A[4] + 593353*A[3] \
                       - 411276*A[2] - 114396*a + 4288) \
            + 48*b*Ap1[1]*(7784*A[6] + 48180*A[5] - 491202*A[4] + 336347*A[3] \
                           + 163734*A[2] - 28908*a - 5560) \
            - (a + 2)*(2*a + 1)*(4568*A[6] - 226668*A[5] - 465702*A[4] \
                                 + 2013479*A[3] - 465702*A[2] - 226668*a \
                                 + 4568)

    C[5] = C[0] / (6688604160.*Ap1[5])
    C[5] *= 1741824*B[10] - 2903040*B[9]*(11*a - 5) \
            + 2177280*B[8]*(110*A[2] - 121*a + 14) \
            - 580608*B[7]*(1628*A[3] - 3333*A[2] + 1023*a + 52) \
            + 169344*B[6]*(12364*A[4] - 43648*A[3] + 26763*A[2] + 1232*a \
                           - 788) \
            - 24192*B[5]*(104852*A[5] - 646624*A[4] + 721391*A[3] \
                          - 16841*A[2] - 74096*a + 148) \
            + 2016*B[4]*(710248*A[6] - 8878716*A[5] + 17928834*A[4] \
                         - 3333407*A[3] - 4339566*A[2] + 287364*a + 89128) \
            - 1344*B[3]*(87824*A[7] - 7150220*A[6] + 29202756*A[5] \
                         - 15113527*A[4] - 14223011*A[3] + 3462492*A[2] \
                         + 1137092*a - 18896) \
            - 84*B[2]*(1690480*A[8] + 14139136*A[7] - 232575464*A[6] \
                       + 296712592*A[5] + 215856619*A[4] - 152181392*A[3] \
                       - 47718440*A[2] + 5813632*a + 943216) \
            + 84*b*Ap1[1]*(82224*A[8] - 5628896*A[7] - 26466520*A[6] \
                           + 168779208*A[5] - 104808005*A[4] - 56259736*A[3] \
                           + 15879912*A[2] + 4020640*a - 63952) \
            + (a + 2)*(2*a + 1)*(2622064*A[8] + 12598624*A[7] \
                                 - 167685080*A[6] - 302008904*A[5] \
                                 + 1115235367*A[4] - 302008904*A[3] \
                                 - 167685080*A[2] + 12598624*a + 2622064)

    Z = pow(a * x, 1/Ap1[1])
    Zp = 1.
    res = C[0]
    for k in range(1, 5):
        Zp /= Z
        res += (-1)**k * C[k] * Zp
    res *= pow(Z, 0.5 - b) * exp(Ap1[1]/a * Z)
    return res


# 4. Integral representation
#    K(a, b, x, r) = exp(-r + x * r^(-a) * cos(pi*a)) * r^(-b)
#                  * sin(x * r^(-a) * sin(pi*a) + pi * b)
#   P(eps, a, b, x, phi) =
#                        * exp(eps * cos(phi) + x * eps^(-a) * cos(a*phi))
#                        * cos(eps * sin(phi) - x * eps^(-a) * sin(a*phi)
#                              + (1-b)*phi)
#
#   Phi(a, b, x) = 1/pi * int_eps^inf K(a, b, x, r) * dr
#        + eps^(1-b)/pi * int_0^pi    P(eps, a, b, x, phi) * dphi
#
#   for any eps > 0.
# Note that P has a misprint in Luchko (2008).
#
# As K has a leading exp(-r), we factor this out and apply Gauss-Laguerre
# quadrature rule:
#   int_0^inf K(a, b, x, r+eps) dr = int_0^inf exp(-r) Kmod(..., r) dr
# The integral over P is done via a Gauss-Legendre quadrature rule.
#
@cython.cdivision(True)
cdef inline double _Kmod(double eps, double a, double b, double x,
                         double r) nogil:
    cdef double x_r_a
    x_r_a = x * pow(r + eps, -a)
    return exp(-eps + x_r_a * cospi(a)) * pow(r + eps, -b) \
        * sin(x_r_a * sinpi(a) + M_PI * b)


@cython.cdivision(True)
cdef inline double _P(double eps, double a, double b, double x,
                      double phi) nogil:
    cdef double eps_a
    x_eps_a = x * pow(eps, -a)
    return exp(eps * cos(phi) + x_eps_a * cos(a*phi)) \
        * cos(eps * sin(phi) - x_eps_a * sin(a*phi) + (1-b)*phi)


@cython.cdivision(True)
cdef inline double wright_bessel_integral(double a, double b, double x) nogil:
    cdef:
        double res1, res2, y, eps
        int k
        # roots of laguerre polynomial of order 50
        # scipy.special.roots_laguerre(50)[0] or
        # sympy.integrals.quadrature.import gauss_laguerre(50, 16)[0]
        double *x_laguerre = \
            [0.02863051833937908, 0.1508829356769337, 0.3709487815348964,
            0.6890906998810479, 1.105625023539913, 1.620961751102501,
            2.23561037591518, 2.950183366641835, 3.765399774405782,
            4.682089387559285, 5.70119757478489, 6.823790909794551,
            8.051063669390792, 9.384345308258407, 10.82510903154915,
            12.37498160875746, 14.03575459982991, 15.80939719784467,
            17.69807093335025, 19.70414653546156, 21.83022330657825,
            24.0791514444115, 26.45405784125298, 28.95837601193738,
            31.59588095662286, 34.37072996309045, 37.28751061055049,
            40.35129757358607, 43.56772026999502, 46.94304399160304,
            50.48426796312992, 54.19924488016862, 58.09682801724853,
            62.18705417568891, 66.48137387844482, 70.99294482661949,
            75.73701154772731, 80.73140480247769, 85.99721113646323,
            91.55969041253388, 97.44956561485056, 103.7048912366923,
            110.3738588076403, 117.5191982031112, 125.2254701334734,
            133.6120279227287, 142.8583254892541, 153.2603719726036,
            165.3856433166825, 180.6983437092145]
        # weights for laguerre polynomial of order 50
        # sympy.integrals.quadrature.import gauss_laguerre(50, 16)[1]
        double *w_laguerre = \
            [0.07140472613518988, 0.1471486069645884, 0.1856716275748313,
            0.1843853825273539, 0.1542011686063556, 0.1116853699022688,
            0.07105288549019586, 0.04002027691150833, 0.02005062308007171,
            0.008960851203646281, 0.00357811241531566, 0.00127761715678905,
            0.0004080302449837189, 0.0001165288322309724, 2.974170493694165e-5,
            6.777842526542028e-6, 1.37747950317136e-6, 2.492886181720092e-7,
            4.010354350427827e-8, 5.723331748141425e-9, 7.229434249182665e-10,
            8.061710142281779e-11, 7.913393099943723e-12, 6.81573661767678e-13,
            5.13242671658949e-14, 3.365624762437814e-15, 1.913476326965035e-16,
            9.385589781827253e-18, 3.950069964503411e-19, 1.417749517827512e-20,
            4.309970276292175e-22, 1.101257519845548e-23, 2.344617755608987e-25,
            4.11854415463823e-27, 5.902246763596448e-29, 6.812008916553065e-31,
            6.237449498812102e-33, 4.452440579683377e-35, 2.426862352250487e-37,
            9.852971481049686e-40, 2.891078872318428e-42, 5.906162708112361e-45,
            8.01287459750397e-48, 6.789575424396417e-51, 3.308173010849252e-54,
            8.250964876440456e-58, 8.848728128298018e-62, 3.064894889844417e-66,
            1.988708229330752e-71, 6.049567152238783e-78]
        # roots of legendre polynomial of order 50
        # sympy.integrals.quadrature.import gauss_legendre(50, 16)[0]
        double *x_legendre = \
            [-0.998866404420071, -0.9940319694320907, -0.9853540840480059,
            -0.9728643851066921, -0.9566109552428079, -0.9366566189448779,
            -0.9130785566557919, -0.885967979523613, -0.8554297694299461,
            -0.8215820708593359, -0.7845558329003993, -0.7444943022260685,
            -0.7015524687068223, -0.6558964656854394, -0.6077029271849502,
            -0.5571583045146501, -0.5044581449074642, -0.4498063349740388,
            -0.3934143118975651, -0.3355002454194374, -0.276288193779532,
            -0.2160072368760418, -0.1548905899981459, -0.09317470156008614,
            -0.03109833832718888, 0.03109833832718888, 0.09317470156008614,
            0.1548905899981459, 0.2160072368760418, 0.276288193779532,
            0.3355002454194374, 0.3934143118975651, 0.4498063349740388,
            0.5044581449074642, 0.5571583045146501, 0.6077029271849502,
            0.6558964656854394, 0.7015524687068223, 0.7444943022260685,
            0.7845558329003993, 0.8215820708593359, 0.8554297694299461,
            0.885967979523613, 0.9130785566557919, 0.9366566189448779,
            0.9566109552428079, 0.9728643851066921, 0.9853540840480059,
            0.9940319694320907, 0.998866404420071]
        # weights for legendre polynomial of order 50
        # sympy.integrals.quadrature.import gauss_legendre(50, 16)[1]
        double *w_legendre = \
            [0.002908622553155141, 0.006759799195745401, 0.01059054838365097,
            0.01438082276148557, 0.01811556071348939, 0.02178024317012479,
            0.02536067357001239, 0.0288429935805352, 0.03221372822357802,
            0.03545983561514615, 0.03856875661258768, 0.0415284630901477,
            0.04432750433880328, 0.04695505130394843, 0.04940093844946632,
            0.05165570306958114, 0.05371062188899625, 0.05555774480621252,
            0.05718992564772838, 0.05860084981322245, 0.05978505870426546,
            0.06073797084177022, 0.06145589959031666, 0.06193606742068324,
            0.06217661665534726, 0.06217661665534726, 0.06193606742068324,
            0.06145589959031666, 0.06073797084177022, 0.05978505870426546,
            0.05860084981322245, 0.05718992564772838, 0.05555774480621252,
            0.05371062188899625, 0.05165570306958114, 0.04940093844946632,
            0.04695505130394843, 0.04432750433880328, 0.0415284630901477,
            0.03856875661258768, 0.03545983561514615, 0.03221372822357802,
            0.0288429935805352, 0.02536067357001239, 0.02178024317012479,
            0.01811556071348939, 0.01438082276148557, 0.01059054838365097,
            0.006759799195745401, 0.002908622553155141]

    # We use the free choice of eps to make the integral better behaved.
    # We try to make eps * sin(phi) and x * eps^(-a) * cos(a*phi) small at
    # the same time. Otherwise, P can be highly oscillatory and unpredictable.
    if x > 1 and a > 1:
        # set eps such that x*eps^(-a) = 1
        eps = pow(x, 1./a)
    elif x > 1 and a > 0:
        # set eps such that eps ~ x * eps^(-rho)
        eps = pow(x, 1./(1+a))
        # safeguard, higher better for larger a, lower better for tiny a.
        eps = min(eps, 100)
    else:
        eps = 1.

    res1 = 0
    res2 = 0
    for k in range(50):
        res1 += w_laguerre[k] * _Kmod(eps, a, b, x, x_laguerre[k])
        # y = (b-a)*(x+1)/2.0 + a  for integration from a=0 to b=pi
        y = M_PI * (x_legendre[k] + 1) / 2.0
        res2 += w_legendre[k] * _P(eps, a, b, x, y)
    # (b-a)/2.0 * np.sum(w*func(y, *args), axis=-1)
    res2 *= M_PI/2.0
    res2 *= pow(eps, 1-b)

    return 1./M_PI * (res1 + res2)


# Note: Hardest argument range is large z, large b and small epsilon
@cython.cdivision(True)
cdef inline double wright_bessel_scalar(double a, double b, double x) nogil:
    cdef:
        double xk_k, res
        int order

    if sc_isnan(a) or sc_isnan(b) or sc_isnan(x):
        return nan
    elif a < 0 or b < 0 or x < 0:
        sf_error.error("wright_bessel", sf_error.DOMAIN, NULL)
        return nan
    elif sc_isinf(x):
        return inf
    elif sc_isinf(a) or sc_isinf(b):
        return 0
    elif b >= rgamma_zero:
          return 0
    elif a >= rgamma_zero:
        # use only first term k=0 of series
        return rgamma(b)
    elif x == 0:
        return rgamma(b)
    elif a == 0:
        return exp(x) * rgamma(b)
    elif (a <= 1e-3 and x <= 1) \
        or (a <= 1e-4 and x <= 10) \
        or (a <= 1e-5 and x <= 100) \
        or (a <= 1e-6 and x < exp_inf):
        # series expansion in 'a' to order=order => precision <~ 1e-13
        # max order = 4
        # a <= 1e-5 and x <= 1    : order = 2
        # a <= 1e-4 and x <= 1    : order = 3
        # a <= 1e-3 and x <= 1    : order = 4
        # a <= 1e-6 and  x <= 10  : order = 2
        # a <= 1e-5 and  x <= 10  : order = 3
        # a <= 1e-6 and  x <= 10  : order = 4
        # a <= 1e-7 and  x <= 100 : order = 2
        # a <= 1e-7 and  x <= 100 : order = 2
        # a <= 1e-6 and  x <= 100 : order = 3
        # a <= 1e-5 and  x <= 100 : order = 4
        order = int(7 + log10(a*x))
        return _wb_small_a(a, b, x, order=order)
    elif x <= 1:
        # 18 term Taylor Series => error mostly smaller 5e-14
        return _wb_series(a, b, x, 0, 18)
    elif x <= 2:
        # 20 term Taylor Series => error mostly smaller 1e-12 to 1e-13
        return _wb_series(a, b, x, 0, 20)
    elif (x >= 500 and a <= 1 and b <= 1) \
        or (x >= 100 and a <= 0.5 and b <= 0.5):
        # asymptotic expansion => precision ~ 8e-11
        return _wb_asymptotic(a, b, x)
    #elif a > 10:
    #    # TODO: Use Taylor Series suitable for large a
    #    return 0
    else:
        return wright_bessel_integral(a, b, x)
