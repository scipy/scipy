# Compute the two-sided one-sample Kolmogorov-Smirnov Prob(Dn <= d) where:
#    D_n = sup_x{|F_n(x) - F(x)|},
#    F_n(x) is the empirical CDF for a sample of size n {x_i: i=1,...,n},
#    F(x) is the CDF of a probability distribution.
#
# Exact methods:
# Prob(D_n >= d) can be computed via a matrix algorithm of Durbin[1]
# or a recursion algorithm due to Pomeranz[2].
# Marsaglia, Tsang & Wang[3] gave a computation-efficient way to perform
# the Durbin algorithm.
# D_n >= d <==>  D_n+ >= d or D_n- >= d (the one-sided K-S statistics), hence
# Prob(D_n >= d) = 2*Prob(D_n+ >= d) - Prob(D_n+ >= d and D_n- >= d).
# For d > 0.5, the latter intersection probability is 0.
#
# Approximate methods:
# For d close to 0.5, ignoring that intersection term may still give a
# reasonable approximation.
# Li-Chien[4] and Korolyuk[5] gave an asymptotic formula extending
# Kolmogorov's initial asymptotic, suitable for large d. (See
# scipy.special.kolmogorov for that asymptotic)
# Pelz-Good[6] used the functional equation for Jacobi theta functions to
# transform the Li-Chien/Korolyuk formula produce a computational formula
# suitable for small d.
#
# Simard and L'Ecuyer[7] provided an algorithm to decide when to use each of
# the above approaches and it is that which is used here.
#
# References:
# [1] Durbin J (1968).
#     "The Probability that the Sample Distribution Function Lies Between Two
#     Parallel Straight Lines."
#     Annals of Mathematical Statistics, 39, 398-411.
# [2] Pomeranz J (1974).
#     "Exact Cumulative Distribution of the Kolmogorov-Smirnov Statistic for
#     Small Samples (Algorithm 487)."
#     Communications of the ACM, 17(12), 703-704.
# [3] Marsaglia G, Tsang WW, Wang J (2003).
#     "Evaluating Kolmogorov's Distribution."
#     Journal of Statistical Software, 8(18), 1-4
# [4] LI-CHIEN, C. (1956).
#     "On the exact distribution of the statistics of A. N. Kolmogorov and
#     their asympt expansion."
#     Acta Matematica Sinica, 6, 55-81.
# [5] KOROLYUK, V. S. (1960).
#     "Asymptotic analysis of the distribution of the maximum deviation in
#     the Bernoulli scheme."
#     Theor. Probability Appl., 4, 339-366.
# [6] Pelz W, Good IJ (1976).
#     "Approximating the Lower Tail-areas of the Kolmogorov-Smirnov One-sample
#     Statistic."
#     Journal of the Royal Statistical Society, Series B, 38(2), 152-156.
#  [7] Simard, R., L'Ecuyer, P. (2011)
# 	"Computing the Two-Sided Kolmogorov-Smirnov Distribution",
# 	Journal of Statistical Software, Vol 39, 11,

import numpy as np
import scipy.special
import scipy.misc

_E128 = 128
_EP128 = np.ldexp(np.longdouble(1), _E128)
_EM128 = np.ldexp(np.longdouble(1), -_E128)

_SQRT2PI = np.sqrt(2 * np.pi)
_MIN_LOG = -708


def _select_and_pin_prob(cdfprob, sfprob, cdf=True):
    """Select either the CDF or SF, and then pin to range 0<=p<=1."""
    p = (cdfprob if cdf else sfprob)
    p = (0.0 if p < 0.0 else (1.0 if p > 1.0 else p))
    return p


def _kolmogn_DMTW(n, d, cdf=True):
    r"""Compute the Kolmogorov CDF:  Pr(D_n <= d) using the MTW approach to
    the Durbin matrix algorithm.

    Marsaglia, Tsang, Wang, 2003."""
    # Write d = (k-h)/n, where k is positive integer and 0 <= h < 1
    # Generate initial matrix H of size m*m where m=(2k-1)
    # Compute k-th row of (n!/n^n) * H^n, scaling intermediate results.
    # Requires memory O(m^2) and computation O(m^2 log(n)).
    # Most suitable for small m.

    if d <= 0.5/n:
        return _select_and_pin_prob(0.0, 1.0, cdf)
    elif d >= 1.0:
        return _select_and_pin_prob(1.0, 0.0, cdf)
    nd = n * d
    k = int(np.ceil(nd))
    h = k - nd
    m = 2 * k - 1

    if k == 1:  # n * d <= 1:
        # p = { 0,                  if h >= 0.5 (0.0 <= d <= 0.5/n)
        #     { (1-2h)^n * n!/n^n,  if h <= 0.5 (0.5/n <= d <= 1/n)
        #        =  n! (2d-1/n)^n
        if n * d <= 0.5:
            return _select_and_pin_prob(0, 1, cdf)
        nu = 2 * d - 1.0 / n
        if n < 10:
            p = np.prod(np.arange(1, n + 1)) * np.power(nu, n)
        else:
            logp = scipy.special.loggamma(n + 1) + n * np.log(nu)
            p = 0.0 if logp < _MIN_LOG else np.exp(logp)
        return _select_and_pin_prob(p, 1.0 - p, cdf)

    # if d >= 0.5, then p = 1 - 2*scipy.special.smirnov(n, d)
    # Pull out the topmost values for d, the ones closest to 1.0
    if k == n:  # d >= 1 - 1/n
        q = 2 * np.power(1 - d, n)  # 2(1-d)^n
        p = 1 - q
        return _select_and_pin_prob(p, q, cdf)

    H = np.zeros([m, m])

    # Initialize: v is first column (and last row) of H
    #  v[j] = (1-h^(j+1)/(j+1)!  (except for v[-1])
    #  w[j] = 1/(j)!
    # q = k-th row of H (actually i!/n^i*H^i)
    intm = np.arange(1, m + 1)
    v = 1.0 - h ** intm
    w = np.zeros(m)
    fac = 1.0
    for j in intm:
        w[j - 1] = fac
        fac /= j  # This might underflow.  Isn't a problem.
        v[j - 1] *= fac
    v[-1] = 1.0 - 2 * h ** m
    if 2 * h - 1.0 > 0:
        v[-1] += (2 * h - 1.0) ** m
    v[-1] *= fac

    for i in range(1, m):
        H[i - 1:, i] = w[:m - i + 1]
    H[:, 0] = v
    H[-1, :] = np.flip(v, axis=0)

    Hpwr = np.eye(np.shape(H)[0])  # Holds intermediate powers of H
    nn = n
    expnt = 0  # Scaling of Hpwr
    Hexpnt = 0  # Scaling of H
    sqr = 1
    while nn > 0:
        if nn % 2:
            Hpwr = np.matmul(Hpwr, H)
            expnt += Hexpnt
        sqr *= 2
        H = np.matmul(H, H)
        Hexpnt *= 2
        # Scale as needed.
        if np.abs(H[k - 1, k - 1]) > _EP128:
            H /= _EP128
            Hexpnt += _E128
        nn = nn // 2

    # ans = (Hpwr[k-1, k-1] * 2^Hpwrexpnt) * n!/n^n
    p = Hpwr[k - 1, k - 1]

    # Multiply by n!/n^n
    for i in range(1, n + 1):
        p = i * p / n
        if np.abs(p) < _EM128:
            p *= _EP128
            expnt -= _E128

    while expnt > 0:
        p = p * _EP128
        expnt -= _E128
    while expnt < 0:
        p = p * _EM128
        expnt += _E128

    assert expnt == 0, str([expnt, p, Hpwr[k - 1, k - 1], Hexpnt])
    return _select_and_pin_prob(p, 1.0-p, cdf)


def _pomeranz_compute_j1j2(i, n, ll, ceilf, roundf):
    '''Compute the endpoints of the interval for row i.'''
    if i == 0:
        j1, j2 = -ll - ceilf - 1, ll + ceilf - 1
    else:
        # i + 1 = 2*idiv2 + remainder
        idiv2 = (i + 1) // 2
        if (i + 1) % 2 == 0:  # i+1 is even, so i is odd
            if idiv2 == n + 1:
                j1, j2 = n - ll - ceilf - 1, n + ll + ceilf - 1
            else:
                j1, j2 = idiv2 - 1 - ll - roundf - 1, idiv2 + ll - 1 + ceilf - 1
        else:
            assert (i + 1) % 2 != 0  # i+1 is odd, i is even
            j1, j2 = idiv2 - 1 - ll - 1, idiv2 + ll + roundf - 1
    return max(j1 + 2, 0), min(j2, n)


def _kolmogn_Pomeranz(n, x, cdf=True):
    r"""Compute the Kolmogorov CDF:  Pr(D_n <= d) using the Pomeranz recursion algorithm."""

    # V is n*(2n+2) matrix.
    # Each row is convolution of previous row and probabilities from a Poisson distribution.
    # Desired prob is n! V[n-1, 2n+2-1]
    # Scale intermediate results as needed
    # Only two rows are needed at any given stage:
    #  - Call them V0 and V1.
    #  - Swap each iteration
    #  - V0s and V1s track the scaling in each row
    # Only a few (contiguous) entries in each row can be non-zero.
    #  - Keep track of start and end (j1 and j2)
    # Only a few different Poisson distributions can occur
    t = n * x
    ll = int(np.floor(t))
    f = t - ll  # fractional part of t
    g = min(f, 1.0 - f)
    ceilf = (1 if f > 0 else 0)
    roundf = (1 if f > 0.5 else 0)
    npwrs = 2 * (ll + 1)    # Maximum number of powers needed in the convolutions
    gpower = np.zeros(npwrs)  # gpower = (g/n)^m/m!
    twogpower = np.zeros(npwrs)  # twogpower = (2g/n)^m/m!
    onem2gpower = np.zeros(npwrs)  # onem2gpower = ((1-2g)/n)^m/m!
    # gpower etc are *almost* Poisson probs, just missing the normalizing factor.
    # This is accounted for

    gpower[0] = 1.0
    twogpower[0] = 1.0
    onem2gpower[0] = 1.0
    expnt = 0
    g_over_n, two_g_over_n, one_minus_two_g_over_n = g/n, 2*g/n, (1 - 2*g)/n
    for m in range(1, npwrs):
        # These may underflow but doesn't matter
        gpower[m] = gpower[m - 1] * g_over_n / m
        twogpower[m] = twogpower[m - 1] * two_g_over_n / m
        onem2gpower[m] = onem2gpower[m - 1] * one_minus_two_g_over_n / m

    V0 = np.zeros([npwrs])
    V1 = np.zeros([npwrs])
    V1[0] = 1  # first row
    V0s, V1s = 0, 0  # scaling of each row

    j1, j2 = _pomeranz_compute_j1j2(0, n, ll, ceilf, roundf)
    for i in range(1, 2 * n + 2):
        # Preserve j1, j2, V1 from last iteration
        k1, k2 = j1, j2
        V0, V1 = V1, V0
        V0s, V1s = V1s, V0s
        V1.fill(0.0)
        j1, j2 = _pomeranz_compute_j1j2(i, n, ll, ceilf, roundf)
        if i == 1 or i == 2 * n + 1:
            pwrs = gpower
        else:
            pwrs = (twogpower if i % 2 else onem2gpower)
        k1m = min(k1, j1)
        assert k1m == k1, str([[n, x], [t, ll, f, g, npwrs], i, [j1, j2], [k1, k2], k1m])
        ln2 = j2 - k1 + 1
        assert ln2 > 0 or (ll == 0 and f <= 0.5), str([[n, x], [t, ll, f, g, npwrs], i, [j1, j2], [k1, k2], ln2])
        if ln2 > 0:
            assert (ln2 <= npwrs), str([[n, x], [t, ll, f, g, npwrs], i, [j1, j2], [k1, k2], ln2])
            assert j2 - j1 + 1 <= npwrs, str([[n, x], [t, ll, f, g, npwrs], i, [j1, j2], [k1, k2], ln2])
            conv = np.convolve(V0[k1 - V0s:k1 - V0s + ln2], pwrs[:ln2])
            assert list(conv[j1 - k1:][:j2 - j1 + 1]) == list(conv[j1 - k1:j1 - k1 + j2 - j1 + 1])
            # V1[:j2 - j1 + 1] = conv[j1 - k1:][:j2 - j1 + 1]
            conv_start = j1 - k1  # First index to use from conv
            conv_len = j2 - j1 + 1  # Number of entries to use from conv
            V1[:conv_len] = conv[conv_start:conv_start + conv_len]
            if 0 < np.max(V1) < _EM128:
                V1 *= _EP128
                expnt -= _E128
            V1s = V0s + j1 - k1
        else:
            V1s = V0s

    # multiply by n!
    ans = V1[n - V1s]
    for m in range(1, n + 1):
        if np.abs(ans) > _EP128:
            ans *= _EM128
            expnt += _E128
        ans *= m

    # Undo any intermediate scaling
    if expnt != 0:
        ans = np.ldexp(ans, expnt)
    ans = _select_and_pin_prob(ans, 1.0 - ans, cdf)
    return ans


def _kolmogn_LCK(n, x, cdf=True):
    """Compute the Li-Chien, Korolyuk approximation to Pr(D_n <= d).

    K0(z) + K1(z)/sqrt(n) + K2(z)/n + K3(z)/n**1.5.

    Li-Chien (1956), Korolyuk (1960)"""
    z = np.sqrt(n)*x
    s = np.zeros(4)
    alpha = -2 * z**2
    zsq = z*z
    zquad = np.power(z, 4)

    # Use a Horner scheme to evaluate sum c_i q^(k^2)
    maxk = int(np.ceil(np.sqrt(_MIN_LOG / alpha)) // 2)
    k = maxk
    while k > 0:
        qpower = np.exp((2 * k - 1) * alpha)
        ksq = k * k
        k4 = k ** 4
        f1 = ksq - (1 if (k % 2) else 0)
        f2 = 5 * ksq + 22 - 15 * (1 if (k % 2) else 0)
        coeffs = np.array(
            [1.0,
             ksq,
             (f1 - 4 * (f1 + 3) * ksq * zsq + 8 * k4 * zquad),
             (f2 / 5 - 4 * (f2 + 45) * ksq * zsq / 15 + 8 * ksq * ksq * zquad) * ksq])
        if k == 0:
            coeffs[0] /= 2.0
        if k % 2:
            coeffs *= -1.0
        s += coeffs
        s *= qpower
        k -= 1

    # Multiply by some constants
    K0to3 = np.array([-2.0, 4 * z / 3.0, 1 / 9.0, -2 * z / 27.0]) * s
    if cdf:
        K0to3 *= -1
        K0to3[0] = 1 - K0to3[0]
    powers_of_n = np.power(n, np.arange(0, len(K0to3)) / 2.0)
    K0to3 = K0to3 / powers_of_n
    Ksum = sum(K0to3)
    return Ksum


def _kolmogn_PelzGood(n, x, cdf=True):
    """Compute the Pelz-Good approximation to Prob(Dn <= x) with 0<=x<=1.

    Start with Li-Chien, Korolyuk approximation:
        K0(z) + K1(z)/sqrt(n) + K2(z)/n + K3(z)/n**1.5
    Transform each K*(z) using Jacobi theta functions into a form suitable
    for small z.

    Pelz-Good 1976."""
    if x <= 0.0:
        return _select_and_pin_prob(0.0, 1.0, cdf=cdf)
    if x >= 1.0:
        return _select_and_pin_prob(1.0, 0.0, cdf=cdf)
    z = np.sqrt(n) * x

    pisquared = np.pi ** 2
    pifour = np.pi ** 4
    pisix = np.pi ** 6

    zsquared = z ** 2
    zfour = z ** 4
    zsix = z ** 6
    q = np.exp(-pisquared / 8 / zsquared)

    # Coefficients of terms in the sums for K1, K2 and K3
    k1a = - zsquared
    k1b = pisquared / 4

    k2a = 6 * zsix + 2 * zfour
    k2b = (2 * zfour - 5 * zsquared) * pisquared / 4
    k2c = pifour * (1 - 2 * zsquared) / 16

    k3d = pisix * (5 - 30 * zsquared) / 64
    k3c = pifour * (-60 * zsquared + 212 * zfour) / 16
    k3b = pisquared * (135 * zfour - 96 * zsix) / 4
    k3a = -30 * zsix - 90 * z ** 8

    K0to3 = np.zeros(4)
    # Sum terms, over odd m
    # Use a Horner scheme to evaluate sum c_i q^(k^2)
    maxk = int(np.ceil(16 * z / np.pi))
    for k in range(maxk, 0, -1):
        m = 2 * k - 1
        msquared = m ** 2
        mfour = m ** 4
        msix = m ** 6
        t = np.power(q, 8 * (k - 1))
        coeffs = np.array([1.0, k1a + k1b * msquared,
                           k2a + k2b*msquared + k2c * mfour,
                           k3a + k3b*msquared + k3c * mfour + k3d * msix])
        K0to3 += coeffs
        K0to3 *= t
    K0to3 *= q

    assert z > 0
    if z > 0:
        K0to3[0] *= _SQRT2PI / z
    else:
        K0to3[0] = 0
    if zfour > 0:
        K0to3[1] *= _SQRT2PI / 6 / zfour
    else:
        K0to3[1] = 0.0
    zseven = z ** 7
    if zseven > 0:
        K0to3[2] *= _SQRT2PI / 72 / zseven
    else:
        K0to3[2] = 0.0
    zten = z ** 10
    if zten > 0:
        K0to3[3] *= _SQRT2PI / 6480 / zten
    else:
        K0to3[3] = 0.0

    # Now do the sum over the other terms, all k
    q = np.exp(-pisquared / 2 / zsquared)
    # Coefficients of terms in the sums for K3
    k3b = 3 * pisquared * zsquared
    k3c = -pifour
    K0to3Extra = np.zeros_like(K0to3)

    for k in range(maxk, 0, -1):
        ksquared = k ** 2
        kfour = k ** 4
        t = np.power(q, 2 * k - 1)
        coeffs = np.array([0.0, 0.0, ksquared, k3b * ksquared + k3c * kfour])
        K0to3Extra += coeffs
        K0to3Extra *= t
    assert z > 0
    zthree = z ** 3
    if zthree > 0:
        K0to3Extra[2] *= -_SQRT2PI * np.pi ** 2 / 36 / zthree
    else:
        K0to3Extra[2] = 0.0
    if zsix > 0:
        K0to3Extra[3] *= _SQRT2PI / 216 / zsix
    else:
        K0to3Extra[3] = 0.0
    K0to3 += K0to3Extra
    powers_of_n = np.power(n * 1.0, np.arange(0, len(K0to3)) / 2.0)
    K0to3 = K0to3 / powers_of_n
    Ksum = sum(K0to3)
    if not cdf:
        Ksum = 1 - Ksum
    return Ksum


def _kolmogn(n, x, cdf=True):
    if int(n) != n or n <= 0:
        return np.nan
    if x >= 1.0:
        return _select_and_pin_prob(1.0, 0.0, cdf=cdf)
    if x <= 0.0:
        return _select_and_pin_prob(0.0, 1.0, cdf=cdf)
    t = n * x
    if t <= 1.0:  # Ruben-Gambino
        if t <= 0.5:
            return _select_and_pin_prob(0.0, 1.0, cdf=cdf)
        prd = 1.0
        mlt = 2 * t - 1
        for m in range(1, n + 1):
            prd = prd * m * mlt / n
        return _select_and_pin_prob(prd, 1.0 - prd, cdf=cdf)
    if t >= n - 1:  # Ruben-Gambino
        onemx = 1.0 - x
        prob = 2 * onemx ** n
        return _select_and_pin_prob(1 - prob, prob, cdf=cdf)
    if x >= 0.5:  # Exact: 2 * smirnov
        prob = 2 * scipy.special.smirnov(n, x)
        return _select_and_pin_prob(1.0 - prob, prob, cdf=cdf)
    nxsquared = t * x
    if n <= 140:
        if nxsquared <= 0.754693:
            prob = _kolmogn_DMTW(n, x, cdf=True)
            return _select_and_pin_prob(prob, 1.0 - prob, cdf=cdf)
        if nxsquared <= 4:
            prob = _kolmogn_Pomeranz(n, x, cdf=True)
            return _select_and_pin_prob(prob, 1.0 - prob, cdf=cdf)
        # Now use Miller approximation of 2*smirnov
        prob = 2 * scipy.special.smirnov(n, x)
        return _select_and_pin_prob(1.0 - prob, prob, cdf=cdf)
    if cdf:
        if n <= 100000:
            if nxsquared >= 18.0:
                return 1.0
            if n * x**1.5 <= 1.4:
                prob = _kolmogn_DMTW(n, x, cdf=True)
                return _select_and_pin_prob(prob, 1.0 - prob, cdf=cdf)
            prob = _kolmogn_PelzGood(n, x, cdf=True)
            return _select_and_pin_prob(prob, 1.0 - prob, cdf=cdf)
        # n > 1e5
        if nxsquared >= 18.0:
            return 1.0
        prob = _kolmogn_PelzGood(n, x, cdf=True)
        return _select_and_pin_prob(prob, 1.0 - prob, cdf=cdf)
    else:
        if nxsquared >= 370.0:
            return 0.0
        if nxsquared >= 2.2:
            prob = 2 * scipy.special.smirnov(n, x)
            return _select_and_pin_prob(1.0 - prob, prob, cdf=cdf)
        cdfprob = _kolmogn(n, x, cdf=True)
        return _select_and_pin_prob(cdfprob, 1.0 - cdfprob, cdf=cdf)


def _kolmogn_p(n, x):
    if int(n) != n or n <= 0:
        return np.nan
    if x >= 1.0 or x <= 0:
        return 0
    t = n * x
    if t <= 1.0:
        if t <= 0.5:
            # Ruben-Gambino
            return 0.0
        prd = 1.0
        mlt = 2 * t - 1
        for m in range(1, n):
            prd = prd * m * mlt / n
        prd = prd * 2 * n * n
        return prd
    if t >= n - 1:
        # Ruben-Gambino
        onemx = 1.0 - x
        pdf = 2 * onemx ** (n-1) * n
        return pdf
    if x >= 0.5:
        pdf = scipy.stats.ksone.pdf(x, n)
        return 2 * pdf

    # Just take a small delta.
    delta = x / 2.0**16
    delta = min(delta, x - 1.0/n)
    delta = min(delta, 0.5 - x)

    def _kk(_x):
        return kolmogn(n, _x)

    pdf = scipy.misc.derivative(_kk, x, dx=delta, order=5)

    return pdf


def kolmogn(n, x, cdf=True):
    """Compute the CDF Pr(Dn <= x) (or its complement), for the two-sided Kolmogorov-Smirnov statistic."""
    it = np.nditer([n, x, None])
    for _n, _x, z in it:
        z[...] = _kolmogn(int(_n), _x, cdf=cdf)
    result = it.operands[2]

    return result


def kolmognp(n, x):
    """Compute the first derivative of Pr(Dn <= x), for the two-sided Kolmogorov-Smirnov statistic."""
    it = np.nditer([n, x, None])
    for _n, _x, z in it:
        z[...] = _kolmogn_p(int(_n), _x)
    result = it.operands[2]

    return result
