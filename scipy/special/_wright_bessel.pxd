# -*-cython-*-
#
# Implementation of Wright's generalized Bessel function Phi, see
# https://dlmf.nist.gov/10.46.E1
#
# Copyright: Christian Lorentzen
#
# Distributed under the same license as SciPy
#
#
# Implementation Overview:
#
# First, different functions are implemented valid for certain domains of the
# three arguments.
# Finally they are put together in wright_bessel_scalar. See the docstring of
# that function for more details.

import cython
from libc.math cimport (cos, exp, floor, fmax, fmin, log, log10, pow, sin,
                        sqrt, isnan, isinf, M_PI, NAN, INFINITY)

cdef extern from "cephes/lanczos.h":
    double lanczos_g

cdef extern from "cephes/polevl.h":
    double polevl(double x, const double coef[], int N) nogil

from ._cephes cimport lgam, rgamma, zeta, sinpi, cospi, lanczos_sum_expg_scaled
from ._digamma cimport digamma
from . cimport sf_error


# rgamma_zero: smallest value x for which rgamma(x) == 0 as x gets large
DEF rgamma_zero = 178.47241115886637
# exp_inf: smallest value x for which exp(x) == inf
DEF exp_inf = 709.78271289338403
# exp_zero: largest value x for which exp(x) == 0
# DEF exp_zero = -745.13321910194117


@cython.cdivision(True)
cdef inline double _exp_rgamma(double x, double y) nogil:
    """Compute exp(x) / gamma(y) = exp(x) * rgamma(y).

    This helper function avoids overflow by using the lanczos approximation
    of the gamma function.
    """
    # gamma = fac * lanczos_sum_expg_scaled(x)
    # fac = ((x + g - 0.5)/e)**(x - 0.5)
    return (exp(x + (1 - log(y + lanczos_g - 0.5)) * (y - 0.5))
            / lanczos_sum_expg_scaled(y))


@cython.cdivision(True)
cdef inline double _wb_series(double a, double b, double x,
    unsigned int nstart, unsigned int nstop) nogil:
    """1. Taylor series expansion in x=0, for x <= 1

    Phi(a, b, x) = sum_k x^k / k! / Gamma(a*k+b)

    Note that every term, and therefore also Phi(a, b, x), is monotone
    decreasing with increasing a or b.
    """
    cdef:
        unsigned int k, k_max
        double xk_k, res

    xk_k = x**nstart * rgamma(nstart + 1)  # x^k/k!
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


@cython.cdivision(True)
cdef inline double _wb_large_a(double a, double b, double x,
    unsigned int n) nogil:
    """2. Taylor series expansion in x=0, for large a.

    Phi(a, b, x) = sum_k x^k / k! / Gamma(a*k+b)

    Use Stirling formula to find k=k_max, the maximum term.
    Then use n terms of Taylor series around k_max.
    """
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


@cython.cdivision(True)
cdef inline double _wb_small_a(double a, double b, double x, int order) nogil:
    """3. Taylor series in a=0 up to order 5, for tiny a and not too large x

    Phi(a, b, x) = exp(x)/Gamma(b)
                   * (1 - a*x * Psi(b) + a^2/2*x*(1+x) * (Psi(b)^2 - Psi'(b)
                      + ... )
                   + O(a^6))

    where Psi is the digamma function.

    Parameter order takes effect only when b > 1e-3 and 2 <= order <= 5,
    otherwise it defaults to 2 or, if b <= 1e-3, to 5. The lower order is, the
    less polygamma functions have to be computed.

    Call: python _precompute/wright_bessel.py 1

    For small b, i.e. b <= 1e-3, cancellation of poles of digamma(b)/Gamma(b)
    and polygamma needs to be carried out => series expansion in a=0 to order 5
    and in b=0 to order 4.
    Call: python _precompute/wright_bessel.py 2
    """
    cdef:
        double dg, pg1, pg2, pg3, pg4, res
        double A[6]  # coefficients of a^k  (1, -x * Psi(b), ...)
        double B[6]  # powers of b^k/k! or terms in polygamma functions
        double C[6]  # coefficients of a^k1 * b^k2
        double X[6]  # polynomials in x
        int k

    if b <= 1e-3:
        # Series expansion of both a and b up to order 5:
        # M_PI = pi
        # M_EG = Euler Gamma aka Euler Mascheroni constant
        # M_Z3 = zeta(3)
        # C[0] = 1
        C[0] = 1.0000000000000000
        # C[1] = 2*M_EG
        C[1] = 1.1544313298030657
        # C[2] = 3*M_EG**2 - M_PI**2/2
        C[2] = -3.9352684291215233
        # C[3] = 4*M_EG**3 - 2*M_EG*M_PI**2 + 8*M_Z3
        C[3] = -1.0080632408182857
        # C[4] = 5*M_EG**4 - 5*M_EG**2*M_PI**2 + 40*M_EG*M_Z3 + M_PI**4/12
        C[4] = 19.984633365874979
        X[0] = 1
        X[1] = x
        X[2] = x*(x + 1)
        X[3] = x*(x*(x + 3) + 1)
        X[4] = x*(x*(x*(x + 6) + 7) + 1)
        X[5] = x*(x*(x*(x*(x + 10) + 25) + 15) + 1)
        B[0] = 1.
        for k in range(1, 5):
            B[k] = b/k * B[k-1]
        # Note that polevl assumes inverse ordering => A[5] = 0th term
        A[5] = rgamma(b)
        A[4] = x         * (C[0] + C[1]*b + C[2]*B[2] + C[3]*B[3] + C[4]*B[4])
        A[3] = X[2]/2.   *        (C[1]   + C[2]*b    + C[3]*B[2] + C[4]*B[3])
        A[2] = X[3]/6.   *                 (C[2]      + C[3]*b    + C[4]*B[2])
        A[1] = X[4]/24.  *                             (C[3]      + C[4]*b)
        A[0] = X[5]/120. *                                          C[4]
        res = exp(x) * polevl(a, A, 5)
        #res = exp(x) * (A[5] + A[4] * a + A[3] * a**2 + A[2] * a**3)
    else:
        # Phi(a, b, x) = exp(x)/gamma(b) * sum(A[i] * X[i] * B[i], i=0..5)
        # A[n] = a^n/n!
        # But here, we repurpose A[n] = X[n] * B[n] / n!
        # Note that polevl assumes inverse ordering => A[order] = 0th term
        dg = digamma(b)
        # pg1 = polygamma(1, b)
        pg1 = zeta(2, b)
        if order <= 2:
            res = 1 + a*x*(-dg + 0.5*a*(1 + x)*(dg**2 - pg1))
        else:
            if order > 5:
                order = 5
            # pg2 = polygamma(2, b)
            pg2 = -2 * zeta(3, b)
            X[0] = 1
            X[1] = x
            X[2] = x*(x + 1)
            X[3] = x*(x*(x + 3) + 1)
            B[0] = 1
            B[1] = -dg
            B[2] = dg**2 - pg1
            B[3] = (-dg**2 + 3*pg1)*dg - pg2
            A[order] = 1
            A[order-1] = X[1] * B[1]
            A[order-2] = X[2] * B[2] / 2.
            A[order-3] = X[3] * B[3] / 6.
            if order >= 4:
                # pg3 = polygamma(3, b)
                pg3 =  6 * zeta(4, b)
                X[4] = x*(x*(x*(x + 6) + 7) + 1)
                B[4] = ((dg**2 - 6*pg1)*dg + 4*pg2)*dg + 3*pg1**2 - pg3
                A[order-4] = X[4] * B[4] / 24.
            if order >= 5:
                # pg4 = polygamma(4, b)
                pg4 = -24 * zeta(5, b)
                X[5] = x*(x*(x*(x*(x + 10) + 25) + 15) + 1)
                B[5] = ((((-dg**2 + 10*pg1)*dg - 10*pg2)*dg - 15*pg1**2
                        + 5*pg3)*dg + 10*pg1*pg2 - pg4)
                A[order-5] = X[5] * B[5] / 120.
            res = polevl(a, A, order)

        # res *= exp(x) * rgamma(b)
        res *= _exp_rgamma(x, b)
    return res


@cython.cdivision(True)
cdef inline double _wb_asymptotic(double a, double b, double x) nogil:
    """4. Asymptotic expansion for large x up to order 8

    Phi(a, b, x) ~ Z^(1/2-b) * exp((1+a)/a * Z) * sum_k (-1)^k * C_k / Z^k

    with Z = (a*x)^(1/(1+a)).
    Call: python _precompute/wright_bessel.py 3
    """
    cdef:
        double A[15]  # powers of a
        double B[17]  # powers of b
        double Ap1[9]  # powers of (1+a)
        double C[9]  # coefficients of asymptotic series a_k
        double Z, Zp, res
        int k

    A[0] = 1.
    B[0] = 1.
    Ap1[0] = 1.
    for k in range(1, 15):
        A[k] = A[k-1] * a
    for k in range(1, 17):
        B[k] = B[k-1] * b
    for k in range(1, 9):
        Ap1[k] = Ap1[k-1] * (1 + a)

    C[0] = 1./sqrt(2. * M_PI * Ap1[1])

    C[1] = C[0] / (24 * Ap1[1])
    C[1] *= (2*a + 1)*(2 + a) - 12*b*(1 + a - b)

    C[2] = C[0] / (1152 * Ap1[2])
    C[2] *= (
        144*B[4]
        - 96*B[3]*(5*a + 1)
        + 24*B[2]*(20*A[2] + 5*a - 4)
        - 24*b*Ap1[1]*(6*A[2] - 7*a - 2)
        + (a + 2)*(2*a + 1)*(2*A[2] - 19*a + 2)
        )

    C[3] = C[0] / (414720 * Ap1[3])
    C[3] *= (
        8640*B[6]
        - 8640*B[5]*(7*a - 1)
        + 10800*B[4]*(14*A[2] - 7*a - 2)
        - 1440*B[3]*(112*A[3] - 147*A[2] - 63*a + 8)
        + 180*B[2]*(364*A[4] - 1288*A[3] - 567*A[2] + 392*a + 76)
        - 180*b*Ap1[1]*(20*A[4] - 516*A[3] + 417*A[2] + 172*a - 12)
        - (a + 2)*(2*a + 1)*(556*A[4] + 1628*A[3] - 9093*A[2] + 1628*a + 556)
        )

    C[4] = C[0] / (39813120 * Ap1[4])
    C[4] *= (
        103680*B[8]
        - 414720*B[7]*(3*a - 1)
        + 725760*B[6]*a*(8*a - 7)
        - 48384*B[5]*(274*A[3] - 489*A[2] + 39*a + 26)
        + 30240*B[4]*(500*A[4] - 1740*A[3] + 495*A[2] + 340*a - 12)
        - 2880*B[3]*(
                     2588*A[5] - 19780*A[4] + 14453*A[3] + 9697*A[2] - 1892*a
                     - 404
                    )
        + 48*B[2]*(
                   11488*A[6] - 547836*A[5] + 1007484*A[4] + 593353*A[3]
                   - 411276*A[2] - 114396*a + 4288
                  )
        + 48*b*Ap1[1]*(
                       7784*A[6] + 48180*A[5] - 491202*A[4] + 336347*A[3]
                       + 163734*A[2] - 28908*a - 5560
                      )
        - (a + 2)*(2*a + 1)*(
                             4568*A[6] - 226668*A[5] - 465702*A[4]
                             + 2013479*A[3] - 465702*A[2] - 226668*a + 4568
                            )
    )

    C[5] = C[0] / (6688604160. * Ap1[5])
    C[5] *= (
        1741824*B[10] - 2903040*B[9]*(11*a - 5)
        + 2177280*B[8]*(110*A[2] - 121*a + 14)
        - 580608*B[7]*(1628*A[3] - 3333*A[2] + 1023*a + 52)
        + 169344*B[6]*(12364*A[4] - 43648*A[3] + 26763*A[2] + 1232*a - 788)
        - 24192*B[5]*(
                      104852*A[5] - 646624*A[4] + 721391*A[3] - 16841*A[2]
                      - 74096*a + 148
                     )
        + 2016*B[4]*(
                     710248*A[6] - 8878716*A[5] + 17928834*A[4] - 3333407*A[3]
                     - 4339566*A[2] + 287364*a + 89128
                    )
        - 1344*B[3]*(
                     87824*A[7] - 7150220*A[6] + 29202756*A[5] - 15113527*A[4]
                     - 14223011*A[3] + 3462492*A[2] + 1137092*a - 18896
                    )
        - 84*B[2]*(
                   1690480*A[8] + 14139136*A[7] - 232575464*A[6]
                   + 296712592*A[5] + 215856619*A[4] - 152181392*A[3]
                   - 47718440*A[2] + 5813632*a + 943216
                  )
        + 84*b*Ap1[1]*(
                       82224*A[8] - 5628896*A[7] - 26466520*A[6]
                       + 168779208*A[5] - 104808005*A[4] - 56259736*A[3]
                       + 15879912*A[2] + 4020640*a - 63952
                      )
        + (a + 2)*(2*a + 1)*(
                             2622064*A[8] + 12598624*A[7]
                             - 167685080*A[6] - 302008904*A[5]
                             + 1115235367.*A[4] - 302008904*A[3]
                             - 167685080*A[2] + 12598624*a + 2622064
                            )
            )

    C[6] = C[0] / (4815794995200. * Ap1[6])
    C[6] *= (
        104509440*B[12] - 209018880*B[11]*(13*a - 7)
        + 574801920*B[10]*(52*A[2] - 65*a + 12)
        - 63866880*B[9]*(2834*A[3] - 6279*A[2] + 2769*a - 134)
        + 23950080*B[8]*(27404*A[4] - 98228*A[3] + 78663*A[2] - 10868*a - 1012)
        - 13685760*B[7]*(
                         105612*A[5] - 599196*A[4] + 791843*A[3] - 224913*A[2]
                         - 27612*a + 4540
                        )
        + 2661120*B[6]*(
                        693680*A[6] - 6473532*A[5] + 13736424*A[4]
                        - 7047469*A[3] - 723840*A[2] + 471588*a + 7376
                       )
        - 2661120*B[5]*(
                        432536*A[7] - 7850804*A[6] + 27531114*A[5]
                        - 24234457*A[4] - 703001*A[3] + 3633474*A[2]
                        - 36244*a - 45128
                       )
        + 166320*B[4]*(
                       548912*A[8] - 75660832*A[7] + 502902712*A[6]
                       - 764807992*A[5] + 91248287*A[4] + 217811464*A[3]
                       - 20365384*A[2] - 9776416*a + 37936
                      )
        + 10080*B[3]*(
                      18759728*A[9] + 165932208*A[8] - 4710418440.*A[7]
                      + 13686052536.*A[6] - 5456818809.*A[5]
                      - 6834514245.*A[4] + 1919299512.*A[3]
                      + 752176152*A[2] - 45661200*a - 8616848
                     )
        - 360*B[2]*(
                    32743360*A[10] - 3381871792.*A[9] - 21488827776.*A[8]
                    + 200389923864.*A[7] - 198708005340.*A[6]
                    - 171633799779.*A[5] + 123124874028.*A[4]
                    + 40072774872.*A[3] - 9137993280.*A[2]
                    - 1895843248.*a + 18929728
                   )
        - 360*b*Ap1[1]*(
                        57685408*A[10] + 406929456*A[9]
                        - 6125375760.*A[8] - 27094918920.*A[7]
                        + 128752249410.*A[6] - 74866710561.*A[5]
                        - 42917416470.*A[4] + 16256951352.*A[3]
                        + 4375268400.*A[2] - 316500688*a - 47197152
                       )
        + (a + 2)*(2*a + 1)*(
                             167898208*A[10] - 22774946512.*A[9]
                             - 88280004528.*A[8] + 611863976472.*A[7]
                             + 1041430242126.*A[6] - 3446851131657.*A[5]
                             + 1041430242126.*A[4] + 611863976472.*A[3]
                             - 88280004528.*A[2] - 22774946512.*a
                             + 167898208
                            )
            )

    C[7] = C[0] / (115579079884800. * Ap1[7])
    C[7] *= (
        179159040*B[14] - 1254113280.*B[13]*(5*a - 3)
        + 1358622720.*B[12]*(70*A[2] - 95*a + 22)
        - 905748480*B[11]*(904*A[3] - 2109*A[2] + 1119*a - 112)
        + 1245404160.*B[10]*(
                             3532*A[4] - 12824*A[3] + 11829*A[2] - 2824*a + 44
                            )
        - 59304960*B[9]*(256820*A[5] - 1397680*A[4] + 2025545*A[3]
                         - 869495*A[2] + 52000*a + 8788)
        + 14826240*B[8]*(2274536*A[6] - 18601572*A[5] + 40698318*A[4]
                         - 28230079*A[3] + 3916398*A[2] + 832668*a - 65176
                        )
        - 59304960*B[7]*(
                         760224*A[7] - 9849164*A[6] + 32495784*A[5]
                         - 34813869*A[4] + 9175207*A[3] + 1898688*A[2]
                         - 469788*a - 13184
                        )
        + 25945920*B[6]*(
                         1167504*A[8] - 28779840*A[7] + 149752856*A[6]
                         - 246026112*A[5] + 111944073*A[4] + 18341600*A[3]
                         - 12131496*A[2] - 274368*a + 102800
                        )
        - 157248*B[5]*(
                       12341872*A[9] - 3122991216.*A[8]
                       + 29900054232.*A[7] - 78024816720.*A[6]
                       + 58914656739.*A[5] + 4637150811.*A[4]
                       - 11523402480.*A[3] + 236218968*A[2] + 337923216*a
                       + 1592048
                      )
        - 28080*B[4]*(
                      265154912*A[10] + 2276098704.*A[9]
                      - 105569461008.*A[8] + 496560666360.*A[7]
                      - 627891462858.*A[6] + 41935358025.*A[5]
                      + 203913875814.*A[4] - 23984801544.*A[3]
                      - 13869306000.*A[2] + 372786832*a + 103532640
                     )
        + 1440*B[3]*(
                     310292864*A[11] - 55169117872.*A[10]
                     - 358957020112.*A[9] + 5714152556088.*A[8]
                     - 13241597459352.*A[7] + 4220720097141.*A[6]
                     + 6845418090249.*A[5] - 2129559215808.*A[4]
                     - 909225098472.*A[3] + 107518582576.*A[2]
                     + 25619444368.*a - 113832704
                    )
        + 12*B[2]*(
                   135319651136.*A[12] + 1119107842176.*A[11]
                   - 22193518174320.*A[10] - 133421793595520.*A[9]
                   + 860103051087996.*A[8] - 703353374803080.*A[7]
                   - 704240127687381.*A[6] + 513111704637960.*A[5]
                   + 166909061348316.*A[4] - 57671564069120.*A[3]
                   - 12453426246000.*A[2] + 695901207936.*a
                   + 93786157376.
                  )
        - 12*b*Ap1[1]*(
                       4365353408.*A[12] - 720248637504.*A[11]
                       - 4222331152560.*A[10] + 29413934270560.*A[9]
                       + 132123980710980.*A[8] - 511247376962820.*A[7]
                       + 283403639131779.*A[6] + 170415792320940.*A[5]
                       - 79274388426588.*A[4] - 21009953050400.*A[3]
                       + 3284035340880.*A[2] + 589294339776.*a
                       - 3693760576.
                      )
        - (a + 2)*(2*a + 1)*(
                             34221025984.*A[12] + 226022948160.*A[11]
                             - 5067505612464.*A[10] - 18868361443936.*A[9]
                             + 86215425028308.*A[8] + 143500920544692.*A[7]
                             - 437682618704613.*A[6] + 143500920544692.*A[5]
                             + 86215425028308.*A[4] - 18868361443936.*A[3]
                             - 5067505612464.*A[2] + 226022948160.*a
                             + 34221025984.
                            )
            )

    C[8] = C[0] / (22191183337881600. * Ap1[8])
    C[8] *= (
        2149908480.*B[16] - 5733089280.*B[15]*(17*a - 11)
        + 7166361600.*B[14]*(272*A[2] - 391*a + 104)
        - 3344302080.*B[13]*(6766*A[3] - 16371*A[2] + 9741*a - 1306)
        + 1811496960.*B[12]*(
                             93092*A[4] - 341564*A[3] + 344199*A[2] - 104924*a
                             + 6308
                            )
        - 517570560*B[11]*(
                           1626220*A[5] - 8641508*A[4] + 13274773*A[3]
                           - 6952303*A[2] + 1007420*a + 5564
                          )
        + 284663808*B[10]*(
                           9979136*A[6] - 75766892*A[5] + 169256148*A[4]
                           - 136824959*A[3] + 35714348*A[2] - 463692*a
                           - 293664
                          )
        - 1423319040.*B[9]*(
                            4466648*A[7] - 49231116*A[6] + 157507414*A[5]
                            - 187114257*A[4] + 78372295*A[3]
                            - 4470082*A[2] - 1913996*a + 82424
                           )
        + 266872320*B[8]*(
                          33133136*A[8] - 564264544*A[7]
                          + 2618606424.*A[6] - 4491310104.*A[5]
                          + 2853943765.*A[4] - 374694552*A[3]
                          - 135365288*A[2] + 17623968*a + 696912
                         )
        - 2156544*B[7]*(
                        2914256144.*A[9] - 93491712432.*A[8]
                        + 664876176984.*A[7] - 1661362937880.*A[6]
                        + 1563719627313.*A[5] - 382840842843.*A[4]
                        - 115399415640.*A[3] + 34565562936.*A[2]
                        + 1609337232.*a - 217321904
                       )
        + 179712*B[6]*(
                       1266018560.*A[10] - 789261834512.*A[9]
                       + 10186841596896.*A[8] - 38877799073352.*A[7]
                       + 54334425968952.*A[6] - 22529574889533.*A[5]
                       - 5132942328000.*A[4] + 3438377465592.*A[3]
                       + 84287641248.*A[2] - 72493479440.*a - 807415936
                      )
        + 13824*B[5]*(
                      156356794976.*A[11] + 1180898077328.*A[10]
                      - 90615270907936.*A[9] + 609258947056248.*A[8]
                      - 1312655191366722.*A[7] + 885900509321745.*A[6]
                      + 112162151855265.*A[5] - 212803071513258.*A[4]
                      + 6805217831352.*A[3] + 10051742651296.*A[2]
                      - 55035924848.*a - 52946379296.
                     )
        - 576*B[4]*(
                    143943926464.*A[12] - 60115486481856.*A[11]
                    - 376366989757200.*A[10] + 9534223075576160.*A[9]
                    - 35603777465262396.*A[8] + 39375990156664980.*A[7]
                    - 868175004137259.*A[6] - 14279180718355020.*A[5]
                    + 1985747535239364.*A[4] + 1264001337603680.*A[3]
                    - 75972792514320.*A[2] - 23855850572736.*a
                    - 4996648256.
                   )
        - 384*B[3]*(
                    2038525473856.*A[13] + 16057322146112.*A[12]
                    - 502133360559024.*A[11] - 2985686417468080.*A[10]
                    + 32418922182093292.*A[9] - 63665380623022452.*A[8]
                    + 16481208821092575.*A[7] + 34161547357596099.*A[6]
                    - 11490298497454932.*A[5] - 5117272758337156.*A[4]
                    + 933703210750480.*A[3] + 234855186762000.*A[2]
                    - 7860524600000.*a - 1226607567040.
                   )
        + 96*B[2]*(
                   324439754752.*A[14] - 77231415197120.*A[13]
                   - 539102931841856.*A[12] + 4618258299956336.*A[11]
                   + 28588485529469792.*A[10] - 141383982651179428.*A[9]
                   + 98783147840417772.*A[8] + 112831723492305801.*A[7]
                   - 83329761150975036.*A[6] - 26553582937192900.*A[5]
                   + 12469117738765952.*A[4] + 2587165396642160.*A[3]
                   - 340406368038080.*A[2] - 53659641606080.*a
                   + 219671272960.
                  )
        + 96*b*Ap1[1]*(
                       1026630779520.*A[14] + 8781958472768.*A[13]
                       - 210659786204384.*A[12] - 1222283505284208.*A[11]
                       + 5064251967491416.*A[10] + 24013052207628140.*A[9]
                       - 79710880160087370.*A[8] + 42596558293213227.*A[7]
                       + 26570293386695790.*A[6] - 14407831324576884.*A[5]
                       - 3617322833922440.*A[4] + 950664948554384.*A[3]
                       + 172358006894496.*A[2] - 7430887938496.*a
                       - 889746675584.
                      )
        - (a + 2)*(2*a + 1)*(
                             573840801152.*A[14]
                             - 156998277198784.*A[13]
                             - 898376974770592.*A[12]
                             + 8622589006459984.*A[11]
                             + 32874204024803560.*A[10]
                             - 111492707520083828.*A[9]
                             - 184768503480287646.*A[8]
                             + 528612016938984183.*A[7]
                             - 184768503480287646.*A[6]
                             - 111492707520083828.*A[5]
                             + 32874204024803560.*A[4]
                             + 8622589006459984.*A[3]
                             - 898376974770592.*A[2]
                             - 156998277198784.*a
                             + 573840801152.
                            )
            )

    Z = pow(a * x, 1/Ap1[1])
    Zp = 1.
    res = C[0]
    for k in range(1, 9):
        Zp /= Z
        res += (-1)**k * C[k] * Zp
    res *= pow(Z, 0.5 - b) * exp(Ap1[1]/a * Z)
    return res


@cython.cdivision(True)
cdef inline double _Kmod(double eps, double a, double b, double x,
                         double r) nogil:
    """Compute integrand Kmod(eps, a, b, x, r) for Gauss-Laguerre quadrature.

    K(a, b, x, r+eps) = exp(-r-eps) * Kmod(eps, a, b, x, r)
    """
    cdef double x_r_a
    x_r_a = x * pow(r + eps, -a)
    return exp(x_r_a * cospi(a)) * pow(r + eps, -b) \
        * sin(x_r_a * sinpi(a) + M_PI * b)


@cython.cdivision(True)
cdef inline double _P(double eps, double a, double b, double x,
                      double phi) nogil:
    """Compute integrand P for Gauss-Legendre quadrature.

    P(eps, a, b, x, phi) =
                         * exp(eps * cos(phi) + x * eps^(-a) * cos(a*phi))
                         * cos(eps * sin(phi) - x * eps^(-a) * sin(a*phi)
                               + (1-b)*phi)
    """
    cdef double eps_a
    x_eps_a = x * pow(eps, -a)
    return exp(eps * cos(phi) + x_eps_a * cos(a*phi)) \
        * cos(eps * sin(phi) - x_eps_a * sin(a*phi) + (1-b)*phi)


@cython.cdivision(True)
cdef inline double wright_bessel_integral(double a, double b, double x) nogil:
    """5. Integral representation

    K(a, b, x, r) = exp(-r + x * r^(-a) * cos(pi*a)) * r^(-b)
                  * sin(x * r^(-a) * sin(pi*a) + pi * b)
    P(eps, a, b, x, phi) =
                         * exp(eps * cos(phi) + x * eps^(-a) * cos(a*phi))
                         * cos(eps * sin(phi) - x * eps^(-a) * sin(a*phi)
                               + (1-b)*phi)

    Phi(a, b, x) = 1/pi * int_eps^inf K(a, b, x, r) * dr
                 + eps^(1-b)/pi * int_0^pi P(eps, a, b, x, phi) * dphi

    for any eps > 0.

    Note that P has a misprint in Luchko (2008).
    This integral representation introduced the free parameter eps (from the
    radius of complex contour integration). We try to choose eps such that
    the integrand behaves smoothly.

    As K has a leading exp(-r), we factor this out and apply Gauss-Laguerre
    quadrature rule:

    int_0^inf K(a, b, x, r+eps) dr = exp(-eps) int_0^inf exp(-r) Kmod(.., r) dr

    Note the shift r -> r+eps to have integation from 0 to infinity.
    The integral over P is done via a Gauss-Legendre quadrature rule.

    Note: Hardest argument range is large z, large b and small eps.
    """
    cdef:
        double res1, res2, y, eps
        int k
        # roots of laguerre polynomial of order 50
        # scipy.special.roots_laguerre(50)[0] or
        # sympy.integrals.quadrature.import gauss_laguerre(50, 16)[0]
        double *x_laguerre = [
            0.02863051833937908, 0.1508829356769337, 0.3709487815348964,
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
            165.3856433166825, 180.6983437092145
            ]
        # weights for laguerre polynomial of order 50
        # sympy.integrals.quadrature.import gauss_laguerre(50, 16)[1]
        double *w_laguerre = [
            0.07140472613518988, 0.1471486069645884, 0.1856716275748313,
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
            1.988708229330752e-71, 6.049567152238783e-78
            ]
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
        double *A = [0.41037, 0.30833, 6.9952, 18.382, -2.8566, 2.1122]

    # We use the free choice of eps to make the integral better behaved.
    # 1. Concern is oscillatory behaviour of P. Therefore, we'd like to
    #    make the change in the argument of cosine small, i.e. make arc length
    #    int_0^phi sqrt(1 + f'(phi)^2) dphi small, with
    #    f(phi) = eps * sin(phi) - x * eps^(-a) * sin(a*phi) + (1-b)*phi
    #    Proxy, make |f'(phi)| small.
    # 2. Concern is int_0 K ~ int_0 (r+eps)^(-b) .. dr
    #    This is difficult as r -> 0  for large b. It behaves better for larger
    #    values of eps.

    # Minimize oscillatory behaviour of P
    eps = (A[0] * b * exp(-0.5 * a)
           + exp(A[1] + 1 / (1 + a) * log(x) - A[2] * exp(-A[3] * a)
                 + A[4] / (1 + exp(A[5] * a))))
    if a >= 4 and x >= 100:
        eps += 1  # This part is hard to fit.

    # Large b
    if b >= 8:
        # Make P small compared to K by setting eps large enough.
        # int K ~ exp(-eps) and int P ~ eps^(1-b)
        eps = max(eps, pow(b, -b / (1. - b)) + 0.1 * b)

    # safeguard, higher better for larger a, lower better for tiny a.
    eps = min(eps, 150.)
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
    """Compute Wright's generalized Bessel function for scalar arguments.

    According to [1], it is an entire function defined as

    .. math:: \Phi(a, b; x) = \sum_{k=0}^\infty \frac{x^k}{k! \Gamma(a k + b)}

    So far, only non-negative values of rho=a, beta=b and z=x are implemented.
    There are 5 different approaches depending on the ranges of the arguments:

    1. Taylor series expansion in x=0 [1], for x <= 1.
       Involves gamma funtions in each term.
    2. Taylor series expansion in x=0 [2], for large a.
    3. Taylor series in a=0, for tiny a and not too large x.
    4. Asymptotic expansion for large x [3, 4].
       Suitable for large x while still small a and b.
    5. Integral representation [5], in principle for all arguments.

    References
    ----------
    [1] https://dlmf.nist.gov/10.46.E1
    [2] P. K. Dunn, G. K. Smyth (2005), Series evaluation of Tweedie exponential
        dispersion model densities. Statistics and Computing 15 (2005): 267-280.
    [3] E. M. Wright (1935), The asymptotic expansion of the generalized Bessel
        function. Proc. London Math. Soc. (2) 38, pp. 257-270.
        https://doi.org/10.1112/plms/s2-38.1.257
    [4] R. B. Paris (2017), The asymptotics of the generalised Bessel function,
        Mathematica Aeterna, Vol. 7, 2017, no. 4, 381 - 406,
        https://arxiv.org/abs/1711.03006
    [5] Y. F. Luchko (2008), Algorithms for Evaluation of the Wright Function for
        the Real Arguments' Values, Fractional Calculus and Applied Analysis 11(1)
        http://sci-gems.math.bas.bg/jspui/bitstream/10525/1298/1/fcaa-vol11-num1-2008-57p-75p.pdf
    """
    cdef:
        double xk_k, res
        int order

    if isnan(a) or isnan(b) or isnan(x):
        return NAN
    elif a < 0 or b < 0 or x < 0:
        sf_error.error("wright_bessel", sf_error.DOMAIN, NULL)
        return NAN
    elif isinf(x):
        if isinf(a) or isinf(b):
            return NAN
        else:
            return INFINITY
    elif isinf(a) or isinf(b):
        return NAN  # or 0
    elif a >= rgamma_zero or b >= rgamma_zero:
        sf_error.error("wright_bessel", sf_error.OVERFLOW, NULL)
        return NAN
    elif x == 0:
        return rgamma(b)
    elif a == 0:
        # return exp(x) * rgamma(b)
        return _exp_rgamma(x, b)
    elif (a <= 1e-3 and b <= 50 and x <= 9) \
        or (a <= 1e-4 and b <= 70 and x <= 100) \
        or (a <= 1e-5 and b <= 170 and x < exp_inf):
        # Taylor Series expansion in a=0 to order=order => precision <= 1e-11
        # If beta is also small => precision <= 1e-11.
        # max order = 5
        if a <= 1e-5:
            if x <= 1:
                order = 2
            elif x <= 10:
                order = 3
            elif x <= 100:
                order = 4
            else:  # x < exp_inf:
                order = 5
        elif a <= 1e-4:
            if x <= 1e-2:
                order = 2
            elif x <= 1:
                order = 3
            elif x <= 10:
                order = 4
            else:  # x <= 100
                order = 5
        else:  # a <= 1e-3
            if x <= 1e-5:
                order = 2
            elif x <= 1e-1:
                order = 3
            elif x <= 1:
                order = 4
            else:  # x <= 9
                order = 5

        return _wb_small_a(a, b, x, order=order)
    elif x <= 1:
        # 18 term Taylor Series => error mostly smaller 5e-14
        return _wb_series(a, b, x, 0, 18)
    elif x <= 2:
        # 20 term Taylor Series => error mostly smaller 1e-12 to 1e-13
        return _wb_series(a, b, x, 0, 20)
    elif a >= 5:
        # Taylor series around the approximate maximum term.
        # Set number of terms=order.
        if a >= 10:
            if x <= 1e11:
                order = 6
            else:
                order = int(fmin(log10(x) - 5 + b/10, 30))
        else:
            if x <= 1e4:
                order = 6
            elif x <= 1e8:
                order = int(2 * log10(x))
            elif x <= 1e10:
                order = int(4 * log10(x) - 16)
            else:
                order = int(fmin(6 * log10(x) - 36, 100))

        return _wb_large_a(a, b, x, order)
    elif (0.5 <= a) & (a <= 1.8) & (100 <= b) & (1e5 <= x):
        return NAN
    elif (pow(a * x, 1 / (1. + a)) >= 14 + b * b / (2 * (1 + a))):
        # Asymptotic expansion in Z = (a*x)^(1/(1+a)) up to 8th term 1/Z^8.
        # For 1/Z^k, the highest term in b is b^(2*k) * a0 / (2^k k! (1+a)^k).
        # As a0 is a common factor to all orders, this explains a bit the
        # domain of good convergence set above.
        # => precision ~ 1e-11 but can go down to ~1e-8 or 1e-7
        # Note: We ensured a <= 5 as this is a bad approximation for large a.
        return _wb_asymptotic(a, b, x)
    else:
        return wright_bessel_integral(a, b, x)
