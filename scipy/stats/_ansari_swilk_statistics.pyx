# cython: boundscheck=False
# cython: initializedcheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: cpow=True

from libc.math cimport exp, sqrt, abs, log, acos
import numpy as np
cimport numpy as cnp
cnp.import_array()


def gscale(int test, int other):
    """
    Cython translation for the FORTRAN 77 code given in:

    Dinneen, L. C. and Blakesley, B. C., "Algorithm AS 93: A Generator for the
    Null Distribution of the Ansari-Bradley W Statistic", Applied Statistics,
    25(1), 1976, :doi:`10.2307/2346534`
    """
    cdef:
        int m = min(test, other)
        int n = max(test, other)
        int astart = ((test + 1) // 2) * (1 + (test // 2))
        int LL = (test * other) // 2 + 1
        int ind
        int n2b1
        int n2b2
        int loop_m = 3
        int part_no
        int ks
        int len1 = 0, len2 = 0, len3 = 0  # Nonzero entry lengths

        bint symm = True if (m+n) % 2 == 0 else False
        bint odd = n % 2
        # We use NumPy arrays for mostly infix operators and slice mechanics
        # Then we define memoryviews on top of these arrays to circumvent
        # the NumPy overhead in the for loop heavy helper functions.
        # In the future, this can be cleaned-up with a better organization
        # of the helper functions. Here we only provided literal FORTRAN
        # translations to remove legacy code.
        cnp.ndarray a1 = cnp.PyArray_ZEROS(1, [LL], cnp.NPY_FLOAT32, 0)
        cnp.ndarray a2 = cnp.PyArray_ZEROS(1, [LL], cnp.NPY_FLOAT32, 0)
        cnp.ndarray a3 = cnp.PyArray_ZEROS(1, [LL], cnp.NPY_FLOAT32, 0)
        float[::1] a1v = a1  # memview of a1
        float[::1] a2v = a2  # memview of a2
        float[::1] a3v = a3  # memview of a3


    if m < 0:
        return 0, np.array([], dtype=np.float32), 2

    # Small cases
    if m == 0:
        a1[0] = 1
        return astart, a1, 0

    if m == 1:
        _start1(a1v, n)
        if not (symm or (other > test)):
            a1v[0], a1v[LL - 1] = 1, 2
        return astart, a1, 0

    if m == 2:
        _start2(a1v, n)
        if not (symm or (other > test)):
            for ind in range(LL//2):
                a1v[LL-1-ind], a1v[ind] = a1v[ind], a1v[LL-1-ind]
        return astart, a1, 0

    with nogil:
        # m > 2, Initialize for odd or even case
        if odd:
            _start1(a1v, n)
            _start2(a2v, n-1)
            len1, len2, n2b1, n2b2 = 1 + (n // 2), n, 1, 2
            part_no = 0
        else:
            _start2(a1v, n)
            _start1(a2v, n-1)
            _start2(a3v, n-2)
            len1, len2, len3, n2b1, n2b2 = (n + 1), n / 2, (n - 1), 2, 1
            part_no = 1


        # loop_m can only increase hence "safe" while
        while loop_m <= m:
            if part_no == 0:
                l1out = _frqadd(a1v, a2v, len2, n2b1)
                len1 += n
                len3 = _imply(a1v, l1out, len1, a3v, loop_m)
                n2b1 += 1
                loop_m += 1
                part_no = 1
            else:
                l2out = _frqadd(a2v, a3v, len3, n2b2)
                len2 += n - 1
                _ = _imply(a2v, l2out, len2, a3v, loop_m)
                n2b2 += 1
                loop_m += 1
                part_no = 0

        if not symm:
            ks = (m + 3) / 2 - 1
            for ind in range(len2):
                a1v[ks+ind] += a2v[ind]

        # reverse array
        if other > test:
            for ind in range(LL//2):
                a1v[LL-1-ind], a1v[ind] = a1v[ind], a1v[LL-1-ind]

    return astart, a1, 0


cdef inline void _start1(float[::1] a, int n) noexcept nogil:
    """
    Helper function for gscale function, see gscale docstring.
    """
    cdef int lout = 1 + (n // 2)

    a[:lout] = 2
    if (n % 2) == 0:
        a[lout-1] = 1


cdef inline void _start2(float[::1] a, int n) noexcept nogil:
    """
    Helper function for gscale function, see gscale docstring.
    """
    cdef:
        int odd = n % 2
        float A = 1.
        float B = 3.
        float C = 2. if odd else 0.
        int ndo = (n + 2 + odd) // 2 - odd

    for ind in range(ndo):
        a[ind] = A
        A += B
        B = 4 - B

    A, B = 1, 3
    for ind in range(n-odd, ndo-1, -1):
        a[ind] = A + C
        A += B
        B = 4 - B

    if odd == 1:
        a[(ndo*2)-1] = 2


cdef inline int _frqadd(float[::1] a, float[::1] b, int lenb,
                        int offset) noexcept nogil:
    """
    Helper function for gscale function, see gscale docstring.
    """
    cdef:
        float two = 2
        int lout = lenb + offset
        int ind

    for ind in range(lenb):
        a[offset+ind] += two*b[ind]
    return lout


cdef int _imply(float[::1] a, int curlen, int reslen, float[::1] b,
                int offset) noexcept nogil:
    """
    Helper function for gscale function, see gscale docstring.
    """
    cdef:
        int i1
        int i2 = -offset
        int j2 = reslen-offset
        int j2min = (j2 + 1) // 2 - 1
        int nextlenb = j2
        int j1 = reslen-1
        float summ
        float diff

    j2 -= 1

    for i1 in range((reslen + 1) // 2):
        if i2 < 0:
            summ = a[i1]
        else:
            summ = a[i1] + b[i2]
            a[i1] = summ

        i2 += 1
        if j2 >= j2min:
            if j1 > curlen - 1:
                diff = summ
            else:
                diff = summ - a[j1]

            b[i1] = diff
            b[j2] = diff
            j2 -= 1

        a[j1] = summ
        j1 -= 1

    return nextlenb


def swilk(double[::1] x, double[::1] a, bint init=False, int n1=-1):
    """
    Calculates the Shapiro-Wilk W test and its significance level

    This is a double precision Cython translation (with modifications) of the
    FORTRAN 77 code given in:

    Royston P., "Remark AS R94: A Remark on Algorithm AS 181: The W-test for
    Normality", 1995, Applied Statistics, Vol. 44, :doi:`10.2307/2986146`

    IFAULT error code details from the R94 paper:
    - 0 for no fault
    - 1 if N, N1 < 3
    - 2 if N > 5000 (a non-fatal error)
    - 3 if N2 < N/2, so insufficient storage for A
    - 4 if N1 > N or (N1 < N and N < 20)
    - 5 if the proportion censored (N-N1)/N > 0.8
    - 6 if the data have zero range (if sorted on input)

    For SciPy, n1 is never used, set to a positive number to enable
    the functionality. Otherwise n1 = n is used.
    """
    cdef:
        int n = len(x)
        int n2 = len(a)
        int ncens, nn2, ind1, i1, ind2
        bint upper = True

        double[6] c1 = [0., 0.221157, -0.147981, -0.207119e1, 0.4434685e1,
                        -0.2706056e1]
        double[6] c2 = [0., 0.42981e-1, -0.293762, -0.1752461e1, 0.5682633e1,
                        -0.3582633e1]
        double[4] c3 = [0.5440, -0.39978, 0.25054e-1, -0.6714e-3]
        double[4] c4 = [0.13822e1, -0.77857, 0.62767e-1, -0.20322e-2]
        double[4] c5 = [-0.15861e1, -0.31082, -0.83751e-1, 0.38915e-2]
        double[3] c6 = [-0.4803, -0.82676e-1, 0.30302e-2]
        double[2] c7 = [0.164, 0.533]
        double[2] c8 = [0.1736, 0.315]
        double[2] c9 = [0.256, -0.635e-2]
        double[2] g = [-0.2273e1, 0.459]
        double Z90 = 0.12816e1
        double Z95 = 0.16449e1
        double Z99 = 0.23263e1
        double ZM = 0.17509e1
        double ZSS = 0.56268
        double BF1 = 0.8378
        double XX90 = 0.556
        double XX95 = 0.622
        double SQRTH = sqrt(2)/2
        double PI6 = 6/np.pi
        double SMALL=1e-19
        double w, pw, an, an25, summ2, ssumm2, rsn
        double A1, A2, fac, delta, w1, y, ld, bf, gamma, m, s
        double RANGE, SA, SX, SSX, SSA, SAX, ASA, XSX, SSASSX, XX, XI
        double Z90F, Z95F, Z99F, ZFM, ZSD, ZBAR

    if n1 < 0:
        n1 = n
    nn2 = n // 2
    if nn2 < n2:
        return 1., 1., 3
    if n < 3:
        return 1., 1., 1
    w = 1.
    pw = 1.
    an = n

    if not init:
        if n == 3:
            a[0] = SQRTH
        else:
            an25 = an + 0.25
            summ2 = 0.
            for ind1 in range(n2):
                temp = _ppnd((ind1+1-0.375)/an25)
                a[ind1] = temp
                summ2 += temp**2

            summ2 *= 2.
            ssumm2 = sqrt(summ2)
            rsn = 1 / sqrt(an)
            A1 = _poly(c1, 6, rsn) - (a[0]/ssumm2)
            if n > 5:
                i1 = 2
                A2 = -a[1]/ssumm2 + _poly(c2, 6, rsn)
                fac = sqrt((summ2 - (2 * a[0]**2) - 2 * a[1]**2) /
                           (1 - (2 * A1**2) - 2 * A2**2))
                a[1] = A2
            else:
                i1 = 1
                fac = sqrt((summ2 - 2 * a[0]**2)/(1 - 2 * A1**2))

            a[0] = A1
            for ind1 in range(i1, nn2):
                a[ind1] *= -1./fac
        init = True

    if n1 < 3:
        return w, pw, 1
    ncens = n - n1

    if ncens < 0 or ((ncens > 0) and (n < 20)):
        return w, pw, 4
    delta = ncens / an
    if delta > 0.8:
        return w, pw, 5

    # Original Fortran code checked the W input value
    # and acted on it via
    #
    # C
    # C If W input as negative, calculate significance level of -W
    # C
    #       IF (W .LT. ZERO) THEN
    #         W1 = ONE + W
    #         IFAULT = 0
    #         GOTO 70
    #       END IF
    #
    # However W was marked as output in the fortran wrappers leading to
    # undefined behavior of the uninitialized W with different compilers.
    # Here W is assumed to be always 0 hence the test is skipped.
    RANGE = x[n1-1] - x[0]
    if RANGE < SMALL:
        return w, pw, 6

    XX = x[0] / RANGE
    SX = XX
    SA = -a[0]
    ind2 = n - 2
    for ind1 in range(1, n1):
        XI = x[ind1] / RANGE
        SX += XI
        if ind1 != ind2:
            SA += (-1 if ind1 < ind2 else 1)*a[min(ind1, ind2)]
        XX = XI
        ind2 -= 1

    ifault = 0
    if n > 5000:
        ifault = 2

    SA /= n1
    SX /= n1
    SSA, SSX, SAX = 0., 0., 0.
    ind2 = n - 1
    for ind1 in range(n1):
        if ind1 != ind2:
            ASA = (-1 if ind1 < ind2 else 1)*a[min(ind1, ind2)] - SA
        else:
            ASA = -SA

        XSX = x[ind1]/RANGE - SX
        SSA += ASA * ASA
        SSX += XSX * XSX
        SAX += ASA * XSX
        ind2 -= 1

    SSASSX = sqrt(SSA * SSX)
    w1 = (SSASSX - SAX) * (SSASSX + SAX)/(SSA * SSX)
    w = 1 - w1

    # Calculate significance level for W (exact for N=3)
    if n == 3:
        # Original Fortran code computation was below
        #
        # pw = PI6 * (asin(sqrt(w)) - PI_OVER_3)
        #
        # However this can return negative p-values for N==3;
        # see gh-18322 and also 32-bit Linux systems.
        # Thus, a potential improvement: precision for small p-values
        # Theoretically w >= 0.75, hence clamping the value
        if w < 0.75:
            return 0.75, 0., ifault
        else:
            pw = 1. - PI6 * acos(sqrt(w))
            return w, pw, ifault

    y = log(w1)
    XX = log(an)
    if n <= 11:
        gamma = _poly(g, 2, an)
        if y >= gamma:
            return w, SMALL, ifault
        y = -log(gamma - y)
        m = _poly(c3, 4, an)
        s = exp(_poly(c4, 4, an))
    else:
        m = _poly(c5, 4, XX)
        s = exp(_poly(c6, 3, XX))

    if ncens > 0:
        ld = -log(delta)
        bf = 1 + XX*BF1
        Z90F = Z90 + bf * _poly(c7, 2, XX90 ** XX)**ld
        Z95F = Z95 + bf * _poly(c8, 2, XX95 ** XX)**ld
        Z99F = Z99 + bf * _poly(c9, 2, XX)**ld
        ZFM = (Z90F + Z95F + Z99F)/3.
        ZSD = (Z90*(Z90F-ZFM)+Z95*(Z95F-ZFM)+Z99*(Z99F-ZFM))/ZSS
        ZBAR = ZFM - ZSD * ZM
        m += ZBAR * s
        s *= ZSD

    pw = _alnorm((y-m)/s, upper)

    return w, pw, ifault


cdef double _alnorm(double x, bint upper) noexcept nogil:
    """
    Helper function for swilk.

    Evaluates the tail area of the standardized normal curve from x to inf
    if upper is True or from -inf to x if upper is False

    Modification has been done to the Fortran version in November 2001 with the
    following note;

        MODIFY UTZERO.  ALTHOUGH NOT NECESSARY
        WHEN USING ALNORM FOR SIMPLY COMPUTING PERCENT POINTS,
        EXTENDING RANGE IS HELPFUL FOR USE WITH FUNCTIONS THAT
        USE ALNORM IN INTERMEDIATE COMPUTATIONS.

    The change is shown below as a commented utzero definition
    """
    cdef:
        double A1, A2, A3, A4, A5, A6, A7
        double B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12
        double ltone = 7.
        # double utzero = 18.66
        double utzero = 38.
        double con = 1.28
        double y, z, temp

    A1 = 0.398942280444
    A2 = 0.399903438504
    A3 = 5.75885480458
    A4 = 29.8213557808
    A5 = 2.62433121679
    A6 = 48.6959930692
    A7 = 5.92885724438
    B1 = 0.398942280385
    B2 = 3.8052e-8
    B3 = 1.00000615302
    B4 = 3.98064794e-4
    B5 = 1.98615381364
    B6 = 0.151679116635
    B7 = 5.29330324926
    B8 = 4.8385912808
    B9 = 15.1508972451
    B10 = 0.742380924027
    B11 = 30.789933034
    B12 = 3.99019417011
    z = x
    if not (z > 0):  # negative of the condition to catch NaNs
        upper = False
        z = -z
    if not ((z <= ltone) or (upper and z <= utzero)):
        return 0. if upper else 1.
    y = 0.5 * z * z
    if z <= con:
        temp = 0.5 - z*(A1 - A2*y/(y + A3 - A4/(y + A5 + A6/(y + A7))))
    else:
        temp = B1 * exp(-y)/(z - B2 + B3/(z + B4 + B5/(z - B6 + B7/
                             (z + B8 - B9/(z + B10 + B11/(z + B12))))))

    return temp if upper else (1-temp)


cdef double _ppnd(double p) noexcept:
    """
    Helper function for swilk. Unused ifault return is
    commented out from return statements.
    """
    cdef:
        double A0, A1, A2, A3, B1, B2, B3, B4, C0, C1, C2, C3, D1, D2
        double q, r, temp
        double split = 0.42

    # cdef int ifault = 0
    A0 = 2.50662823884
    A1 = -18.61500062529
    A2 = 41.39119773534
    A3 = -25.44106049637
    B1 = -8.47351093090
    B2 = 23.08336743743
    B3 = -21.06224101826
    B4 = 3.13082909833
    C0 = -2.78718931138
    C1 = -2.29796479134
    C2 = 4.85014127135
    C3 = 2.32121276858
    D1 = 3.54388924762
    D2 = 1.63706781897

    q = p - 0.5
    if abs(q) <= split:
        r = q * q
        temp = q*(((A3*r + A2)*r + A1) * r + A0)
        temp = temp / ((((B4*r + B3)*r + B2)*r + B1)*r + 1.)
        return temp  #, 0

    r = p
    if q > 0:
        r = 1 - p
    if r > 0:
        r = sqrt(-log(r))
    else:
        return 0.  #, 1

    temp = (((C3*r + C2)*r + C1)*r + C0)
    temp /= (D2*r + D1)*r + 1.
    return (-temp if q < 0 else temp)  #, 0


cdef double _poly(double[::1]c, int nord, double x) noexcept nogil:
    """
    Helper function for swilk function that evaluates polynomials.
    For some reason, the coefficients are given as

        [c0, cn, cn-1, ..., c2, c1]

    hence the backwards loop.
    """
    cdef double res, p
    cdef int ind
    res = c[0]
    if nord == 1:
        return res

    p = x*c[nord-1]
    if nord == 2:
        return res + p

    for ind in range(nord-2, 0, -1):
        p = (p + c[ind])*x
    res += p
    return res
