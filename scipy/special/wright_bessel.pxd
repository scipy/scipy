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
# There are 5 different approaches depending on the ranges of the arguments.
# 1. Taylor series expansion in x=0 [1], for x <= 1.
#    Involves gamma funtions in each term.
# 2. Taylor series expansion in x=0 [2], for large a.
# 3. Taylor series in a=0, for tiny a and not too large x.
# 4. Asymptotic expansion for large x [3, 4].
#    Suitable for large x while still small a and b.
# 5. Integral representation [5], in principle for all arguments.
#
# References:
# [1] https://dlmf.nist.gov/10.46.E1
# [2] P. K. Dunn, G. K. Smyth (2005), Series evaluation of Tweedie exponential
#     dispersion model densities. Statistics and Computing 15 (2005): 267-280.
# [3] E. M. Wright (1935), The asymptotic expansion of the generalized Bessel
#     function. Proc. London Math. Soc. (2) 38, pp. 257–270.
#     https://doi.org/10.1112/plms/s2-38.1.257
# [4] R. B. Paris (2017), The asymptotics of the generalised Bessel function,
#     Mathematica Aeterna, Vol. 7, 2017, no. 4, 381 - 406,
#     https://arxiv.org/abs/1711.03006
# [5] Y. F. Luchko (2008), Algorithms for Evaluation of the Wright Function for
#     the Real Arguments’ Values, Fractional Calculus and Applied Analysis 11(1)
#     http://sci-gems.math.bas.bg/jspui/bitstream/10525/1298/1/fcaa-vol11-num1-2008-57p-75p.pdf

import cython
from libc.math cimport cos, exp, floor, fmax, log, log10, pow, sin, sqrt, M_PI

cdef extern from "_c99compat.h":
    int sc_isnan(double x) nogil
    int sc_isinf(double x) nogil

from ._cephes cimport Gamma, lgam, rgamma, zeta, sinpi, cospi
from ._digamma cimport digamma
from ._complexstuff cimport inf, nan
from . cimport sf_error

# Euler Gamma constant aka Euler Mascheroni constant
DEF M_EG = 0.57721566490153286060
# zeta(3)
DEF M_Z3 = 1.2020569031595942854
# rgamma_zero: smallest value x for which rgamma(x) == 0 as x gets large
DEF rgamma_zero = 178.47241115886637
# exp_inf: smallest value x for which exp(x) == inf
DEF exp_inf = 709.78271289338403
# exp_zero: largest value x for which exp(x) == 0
# DEF exp_zero = -745.13321910194117


# 1. Taylor series expansion in x=0, for x <= 1
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


# 2. Taylor series expansion in x=0, for large a.
#    Phi(a, b, x) = sum_k x^k / k! / Gamma(a*k+b)
#
#    Use Stirling formula to find k=k_max, the maximum term.
#    Then use n terms of Taylor series around k_max.
@cython.cdivision(True)
cdef inline double _wb_large_a(double a, double b, double x,
    unsigned int n) nogil:
    cdef:
        unsigned int k, k_max
        int n_start
        double lnx, res

    k_max = int(pow(pow(a, -a) * x, 1./(1 + a)))
    n_start = k_max - n//2
    if n_start < 0:
        n_start = 0

    res = 0
    lnx = log(x)
    for k in range(n_start, n_start + n):
        res += exp(k * lnx - lgam(k + 1) - lgam(a*k + b))

    return res


# 3. Taylor series in a=0 up to order 5, for tiny a and not too large x
# Phi(a, b, x) = exp(x)/Gamma(b)
#               * (1 - a*x*Psi(b) + a^2/2*x*(1+x)(Psi(b)^2 - Psi'(b)) + O(a^3))
# where Psi is the digamma function.
# Call: python _precompute/wright_bessel.py 1
#
# For small b, cancellation of poles of digamma(b)/Gamma(b) and polygamma needs
# to be carried out => series expansion in a and b around 0 to order 5.
# Call: python _precompute/wright_bessel.py 2
@cython.cdivision(True)
cdef inline double _wb_small_a(double a, double b, double x, int order) nogil:
    cdef:
        double dg, dg1, dg2, dg3, res
        double A[6]  # powers of a^k/k!
        double B[5]  # powers of b^k/k!
        double C[5]  # coefficients of a^k1 * b^k2
        int k

    if b <= 1e-3:
        # Series expansion of both a and b up to order 5:
        C[0] = 1
        C[1] = 2*M_EG
        C[2] = 3*M_EG**2 - M_PI**2/2
        C[3] = 4*M_EG**3 - 2*M_EG*M_PI**2 + 8*M_Z3
        C[4] = 5*M_EG**4 - 5*M_EG**2*M_PI**2 + 40*M_EG*M_Z3 + M_PI**4/12
        A[0] = 1.
        B[0] = 1.
        for k in range(1, 5):
            A[k] = a/k * A[k-1]
            B[k] = b/k * B[k-1]
        A[5] = a/5*A[4]
        res = rgamma(b)
        res += a*x * (C[0] + C[1]*b + C[2]*B[2] + C[3]*B[3] + C[4]*B[4])
        res += A[2]*x*(1+x)*(C[1]   + C[2]*b    + C[3]*B[2] + C[4]*B[3])
        res += A[3]*x*(x**2+3*x+1) * (C[2]      + C[3]*b    + C[4]*B[2])
        res += A[4]*x*(x**3+6*x**2+7*x+1) *      (C[3]      + C[4]*b)
        res += A[5]*x*(x**4+10*x**3+25*x**2+15*x+1) *         C[4]
        res *= exp(x)
    else:
        dg = digamma(b)
        # dg1 = polygamma(1, b)
        dg1 = zeta(2, b)
        # res = 1 - a*x*dg + a**2/2*x*(1+x)*(dg**2 - dg1)
        res = 1 + a*x*(-dg + 0.5*a*(1 + x)*(dg**2 - dg1))
        if order >= 3:
            # dg2 = polygamma(2, b)
            dg2 = -2. * zeta(3, b)
            res += -a**3/6. * x*(x**2 + 3*x + 1)*(dg**3 - 3*dg*dg1 + dg2)
        if order >= 4:
            # dg3 = polygamma(3, b)
            dg3 =  6 * zeta(4, b)
            res += a**4/24. * x*(x**3 + 6*x**2 + 7*x + 1) \
                   *(dg**4 - 6*dg**2*dg1 + 4*dg*dg2 + 3*dg1**2 - dg3)
        res *= exp(x) * rgamma(b)
    return res


# 4. Asymptotic expansion for large x up to order 6
#    Phi(a, b, x) ~ Z^(1/2-b) * exp((1+a)/a * Z) * sum_k (-1)^k * C_k / Z^k
#    with Z = (a*x)^(1/(1+a)).
#
# Call: python _precompute/wright_bessel.py 3
@cython.cdivision(True)
cdef inline double _wb_asymptotic(double a, double b, double x) nogil:
    cdef:
        double A[11]  # powers of a
        double B[13]  # powers of b
        double Ap1[7]  # powers of (1+a)
        double C[7]  # coefficients of asymptotic series a_k
        double Z, Zp, res
        int k

    A[0] = 1.
    B[0] = 1.
    Ap1[0] = 1.
    for k in range(1, 11):
        A[k] = A[k-1] * a
    for k in range(1, 13):
        B[k] = B[k-1] * b
    for k in range(1, 7):
        Ap1[k] = Ap1[k-1] * (1 + a)

    C[0] = 1./sqrt(2. * M_PI * Ap1[1])

    C[1] = C[0] / (24 * Ap1[1])
    C[1] *= (2*a + 1)*(2 + a) - 12*b*(1 + a - b)

    C[2] = C[0] / (1152 * Ap1[2])
    C[2] *= 144*B[4] - 96*B[3]*(5*a + 1) + 24*B[2]*(20*A[2] + 5*a - 4) \
            - 24*b*Ap1[1]*(6*A[2] - 7*a - 2) \
            + (a + 2)*(2*a + 1)*(2*A[2] - 19*a + 2)

    C[3] = C[0] / (414720 * Ap1[3])
    C[3] *= 8640*B[6] - 8640*B[5]*(7*a - 1) + 10800*B[4]*(14*A[2] - 7*a - 2) \
            - 1440*B[3]*(112*A[3] - 147*A[2] - 63*a + 8) \
            + 180*B[2]*(364*A[4] - 1288*A[3] - 567*A[2] + 392*a + 76) \
            - 180*b*Ap1[1]*(20*A[4] - 516*A[3] + 417*A[2] + 172*a - 12) \
            - (a + 2)*(2*a + 1)*(556*A[4] + 1628*A[3] - 9093*A[2] + 1628*a \
                                 + 556)

    C[4] = C[0] / (39813120 * Ap1[4])
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
                       + 2013479*A[3] - 465702*A[2] - 226668*a + 4568)

    C[5] = C[0] / (6688604160. * Ap1[5])
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
                                 + 1115235367.*A[4] - 302008904*A[3] \
                                 - 167685080*A[2] + 12598624*a + 2622064)

    C[6] = C[0] / (4815794995200. * Ap1[6])
    C[6] *= 104509440*B[12] - 209018880*B[11]*(13*a - 7) \
            + 574801920*B[10]*(52*A[2] - 65*a + 12) \
            - 63866880*B[9]*(2834*A[3] - 6279*A[2] + 2769*a - 134) \
            + 23950080*B[8]*(27404*A[4] - 98228*A[3] + 78663*A[2] - 10868*a \
                             - 1012) \
            - 13685760*B[7]*(105612*A[5] - 599196*A[4] + 791843*A[3] \
                             - 224913*A[2] - 27612*a + 4540) \
            + 2661120*B[6]*(693680*A[6] - 6473532*A[5] + 13736424*A[4] \
                            - 7047469*A[3] - 723840*A[2] + 471588*a + 7376) \
            - 2661120*B[5]*(432536*A[7] - 7850804*A[6] + 27531114*A[5] \
                            - 24234457*A[4] - 703001*A[3] + 3633474*A[2] \
                            - 36244*a - 45128) \
            + 166320*B[4]*(548912*A[8] - 75660832*A[7] + 502902712*A[6] \
                           - 764807992*A[5] + 91248287*A[4] + 217811464*A[3] \
                           - 20365384*A[2] - 9776416*a + 37936) \
            + 10080*B[3]*(18759728*A[9] + 165932208*A[8] - 4710418440.*A[7] \
                          + 13686052536.*A[6] - 5456818809.*A[5] \
                          - 6834514245.*A[4] + 1919299512.*A[3] \
                          + 752176152*A[2] - 45661200*a - 8616848) \
            - 360*B[2]*(32743360*A[10] - 3381871792.*A[9] - 21488827776.*A[8] \
                        + 200389923864.*A[7] - 198708005340.*A[6] \
                        - 171633799779.*A[5] + 123124874028.*A[4] \
                        + 40072774872.*A[3] - 9137993280.*A[2] \
                        - 1895843248.*a + 18929728) \
            - 360*b*Ap1[1]*(57685408*A[10] + 406929456*A[9] \
                            - 6125375760.*A[8] - 27094918920.*A[7] \
                            + 128752249410.*A[6] - 74866710561.*A[5] \
                            - 42917416470.*A[4] + 16256951352.*A[3] \
                            + 4375268400.*A[2] - 316500688*a - 47197152) \
            + (a + 2)*(2*a + 1)*(167898208*A[10] - 22774946512.*A[9] \
                                 - 88280004528.*A[8] + 611863976472.*A[7] \
                                 + 1041430242126.*A[6] - 3446851131657.*A[5] \
                                 + 1041430242126.*A[4] + 611863976472.*A[3] \
                                 - 88280004528.*A[2] - 22774946512.*a \
                                 + 167898208)

    Z = pow(a * x, 1/Ap1[1])
    Zp = 1.
    res = C[0]
    for k in range(1, 7):
        Zp /= Z
        res += (-1)**k * C[k] * Zp
    res *= pow(Z, 0.5 - b) * exp(Ap1[1]/a * Z)
    return res


# 5. Integral representation
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
#   int_0^inf K(a, b, x, r+eps) dr = exp(-eps) int_0^inf exp(-r) Kmod(.., r) dr
# Note the shift r -> r+eps to have integation from 0 to infinity.
#
# The integral over P is done via a Gauss-Legendre quadrature rule.
@cython.cdivision(True)
cdef inline double _Kmod(double eps, double a, double b, double x,
                         double r) nogil:
    cdef double x_r_a
    x_r_a = x * pow(r + eps, -a)
    return exp(x_r_a * cospi(a)) * pow(r + eps, -b) \
        * sin(x_r_a * sinpi(a) + M_PI * b)


@cython.cdivision(True)
cdef inline double _P(double eps, double a, double b, double x,
                      double phi) nogil:
    cdef double eps_a
    x_eps_a = x * pow(eps, -a)
    return exp(eps * cos(phi) + x_eps_a * cos(a*phi)) \
        * cos(eps * sin(phi) - x_eps_a * sin(a*phi) + (1-b)*phi)


# Note: Hardest argument range is large z, large b and small epsilon.
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
        # Fitted parameters for optimal choice of eps
        # Call: python _precompute/wright_bessel.py 4
        double *A = [0.26575, 55.827, 0.62564, 1.2511, 4.7228, 8.8881e-05]

    # We use the free choice of eps to make the integral better behaved.
    # 1. Concern is int_0 K ~ int_0 (r+eps)^(-b) .. dr
    #    This is difficult as r -> 0  for large b. It behaves better for larger
    #    values of eps.
    # 2. Concern is oscillatory behaviour of P. Therefore, we'd like to
    #    make the change in the argument of cosine small, i.e. make arc length
    #    int_0^phi sqrt(1 + f'(phi)^2) dphi small, with
    #    f(phi) = eps * sin(phi) - x * eps^(-a) * sin(a*phi) + (1-b)*phi
    #    Proxy, make |f'(phi)| small.
    if b >= 8:
        # Make P small compared to K by setting eps large enough.
        # int K ~ exp(-eps) and int P ~ eps^(1-b)
        eps = pow(b, -b / (1. - b))
    else:
        eps = (A[0] * b + A[1] * pow(x, A[2]) * pow(a, A[3])
               * (exp(-A[4] * sqrt(a)) + A[5]))
    # safeguard, higher better for larger a, lower better for tiny a.
    eps = min(eps, 50.)
    eps = max(eps, 3.)  # Note: 3 seems to be a pretty good choice in general.

    res1 = 0
    res2 = 0
    for k in range(50):
        res1 += w_laguerre[k] * _Kmod(eps, a, b, x, x_laguerre[k])
        # y = (b-a)*(x+1)/2.0 + a  for integration from a=0 to b=pi
        y = M_PI * (x_legendre[k] + 1) / 2.0
        res2 += w_legendre[k] * _P(eps, a, b, x, y)
    res1 *= exp(-eps)
    # (b-a)/2.0 * np.sum(w*func(y, *args), axis=-1)
    res2 *= M_PI/2.0
    res2 *= pow(eps, 1-b)

    return 1./M_PI * (res1 + res2)


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
        if sc_isinf(a) or sc_isinf(b):
            return nan
        else:
            return inf
    elif sc_isinf(a) or sc_isinf(b):
        return nan  # or 0
    elif b >= rgamma_zero:
        sf_error.error("wright_bessel", sf_error.OVERFLOW, NULL)
        return nan # or 0
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
        # if beta also small => precision <~ 2e-14
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
    elif (a >= 10) or (a >= 5 and x <= 1e9):
        # 20 term Taylor series around the approximate maximum term
        # => error mostly smaller 1e-13
        return _wb_large_a(a, b, x, 20)
    elif (a <= 2) and (a * x >= pow(12. * fmax(1., b), 1 + a)):
        # asymptotic expansion in Z = (a*x)^(1/(1+a)) >= 12 * max(1, b)
        # => precision ~ 1e-11 but can go down to ~1e-8 or 1e-7
        return _wb_asymptotic(a, b, x)
    else:
        return wright_bessel_integral(a, b, x)
