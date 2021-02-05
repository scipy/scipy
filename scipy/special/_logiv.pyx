import cython

from libc.math cimport sqrt, log, log1p, M_PI


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef inline double log_iv_asym(double v, double z) nogil:
    r"""Asymptotic expansion for I_{v}(z)
    for real $z > 0$ and $v\to +\infty$.
    Based off DLMF 10.41

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/10.41#ii
    .. [a144617] Sequence A144617, The On-Line Encyclopedia of Integer Sequences,
           https://oeis.org/A144617
    """
    cdef:
        # DLMF 10.41.3 calculates I_v(v*z)
        double x = z / v
        # DLMF 10.41.8
        double inv_p = sqrt(1.0 + x*x)
        # DLMF 10.41.7
        double eta = inv_p + log(x) - log1p(inv_p)
        # e^(vn) / ((2piv)^(1/2)*p^(-1/2))
        double log_res = v*eta - log(2.0 * M_PI * v * inv_p) / 2

        double p, p2, p4, p6, p8, p10
        double pv, pv2, pv3, pv4, pv5

    # large-v asymptotic correction, DLMF 10.41.10
    p = 1.0 / inv_p
    p2 = p * p
    p4 = p2 * p2
    p6 = p4 * p2
    p8 = p4 * p4
    p10 = p6 * p4

    # p over v
    pv  = p / v
    pv2 = pv * pv
    pv3 = pv2 * pv
    pv4 = pv2 * pv2
    pv5 = pv3 * pv2

    # terms from OEIS A144617 and A144618
    u1 = (3.0 - 5.0*p2) / 24.0
    u2 = (81.0 - 462.0*p2 + 385.0*p4) / 1152.0
    u3 = (30375 - 369603*p2 + 765765*p4 - 425425*p6) / 414720.0
    u4 = (4465125 - 94121676*p2 + 349922430*p4 - 446185740*p6 +
          185910725*p8) / 39813120.0
    u5 = (1519035525.0 - 49286948607.0*p2 + 284499769554.0*p4 - 614135872350.0*p6 +
          566098157625.0*p8 - 188699385875.0*p10) / 6688604160.0
    u_corr_i = 1.0 + u1 * pv + u2 * pv2 + u3 * pv3 + u4 * pv4 + u5 * pv5

    return log_res + log(u_corr_i)
