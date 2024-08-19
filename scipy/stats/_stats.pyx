# cython: cpow=True

from cpython cimport bool
from libc cimport math
from libc.math cimport NAN, INFINITY, M_PI as PI
cimport cython
cimport numpy as np
from numpy cimport ndarray, int64_t, float64_t, intp_t

import warnings
import numpy as np
import scipy.stats, scipy.special
from scipy.linalg import solve_triangular
cimport scipy.special.cython_special as cs

np.import_array()


cdef double von_mises_cdf_series(double k, double x, unsigned int p) noexcept:
    cdef double s, c, sn, cn, R, V
    cdef unsigned int n
    s = math.sin(x)
    c = math.cos(x)
    sn = math.sin(p * x)
    cn = math.cos(p * x)
    R = 0
    V = 0
    for n in range(p - 1, 0, -1):
        sn, cn = sn * c - cn * s, cn * c + sn * s
        R = k / (2 * n + k * R)
        V = R * (sn / n + V)

    with cython.cdivision(True):
        return 0.5 + x / (2 * PI) + V / PI


cdef von_mises_cdf_normalapprox(k, x):
    cdef double SQRT2_PI = 0.79788456080286535588  # sqrt(2/pi)

    b = SQRT2_PI / scipy.special.i0e(k)  # Check for negative k
    z = b * np.sin(x / 2.)
    return scipy.stats.norm.cdf(z)


@cython.boundscheck(False)
def von_mises_cdf(k_obj, x_obj):
    cdef double[:] temp, temp_xs, temp_ks
    cdef unsigned int i, p
    cdef double a1, a2, a3, a4, CK
    cdef np.ndarray k = np.asarray(k_obj)
    cdef np.ndarray x = np.asarray(x_obj)
    cdef bint zerodim = k.ndim == 0 and x.ndim == 0

    k = np.atleast_1d(k)
    x = np.atleast_1d(x)
    ix = np.round(x / (2 * PI))
    x = x - ix * (2 * PI)

    # These values should give 12 decimal digits
    CK = 50
    a1, a2, a3, a4 = 28., 0.5, 100., 5.

    bx, bk = np.broadcast_arrays(x, k)
    result = np.empty_like(bx, float)

    c_small_k = bk < CK
    temp = result[c_small_k]
    temp_xs = bx[c_small_k].astype(float)
    temp_ks = bk[c_small_k].astype(float)
    for i in range(len(temp)):
        p = <int>(1 + a1 + a2 * temp_ks[i] - a3 / (temp_ks[i] + a4))
        temp[i] = von_mises_cdf_series(temp_ks[i], temp_xs[i], p)
        temp[i] = 0 if temp[i] < 0 else 1 if temp[i] > 1 else temp[i]
    result[c_small_k] = temp
    result[~c_small_k] = von_mises_cdf_normalapprox(bk[~c_small_k], bx[~c_small_k])

    if not zerodim:
        return result + ix
    else:
        return (result + ix)[0]

@cython.wraparound(False)
@cython.boundscheck(False)
def _kendall_dis(intp_t[:] x, intp_t[:] y):
    cdef:
        intp_t sup = 1 + np.max(y)
        # Use of `>> 14` improves cache performance of the Fenwick tree (see gh-10108)
        intp_t[::1] arr = np.zeros(sup + ((sup - 1) >> 14), dtype=np.intp)
        intp_t i = 0, k = 0, size = x.size, idx
        int64_t dis = 0

    with nogil:
        while i < size:
            while k < size and x[i] == x[k]:
                dis += i
                idx = y[k]
                while idx != 0:
                    dis -= arr[idx + (idx >> 14)]
                    idx = idx & (idx - 1)

                k += 1

            while i < k:
                idx = y[i]
                while idx < sup:
                    arr[idx + (idx >> 14)] += 1
                    idx += idx & -idx
                i += 1

    return dis


# The weighted tau will be computed directly between these types.
# Arrays of other types will be turned into a rank array using _toint64().

ctypedef fused ordered:
    np.int32_t
    np.int64_t
    np.float32_t
    np.float64_t


# Inverts a permutation in place [B. H. Boonstra, Comm. ACM 8(2):104, 1965].
@cython.wraparound(False)
@cython.boundscheck(False)
cdef _invert_in_place(intp_t[:] perm):
    cdef intp_t n, i, j, k
    for n in range(len(perm)-1, -1, -1):
        i = perm[n]
        if i < 0:
            perm[n] = -i - 1
        else:
            if i != n:
                k = n
                while True:
                    j = perm[i]
                    perm[i] = -k - 1
                    if j == n:
                        perm[n] = i
                        break

                    k = i
                    i = j


@cython.wraparound(False)
@cython.boundscheck(False)
def _toint64(x):
    cdef intp_t i = 0, j = 0, l = len(x)
    cdef intp_t[::1] perm = np.argsort(x, kind='quicksort')
    # The type of this array must be one of the supported types
    cdef int64_t[::1] result = np.ndarray(l, dtype=np.int64)

    # Find nans, if any, and assign them the lowest value
    for i in range(l - 1, -1, -1):
        if not np.isnan(x[perm[i]]):
            break
        result[perm[i]] = 0

    if i < l - 1:
        j = 1
        l = i + 1

    for i in range(l - 1):
        result[perm[i]] = j
        if x[perm[i]] != x[perm[i + 1]]:
            j += 1

    result[perm[l - 1]] = j
    return np.array(result, dtype=np.int64)


@cython.wraparound(False)
@cython.boundscheck(False)
def _weightedrankedtau(const ordered[:] x, const ordered[:] y, intp_t[:] rank, weigher, bool additive):
    # y_local and rank_local (declared below) are a work-around for a Cython
    # bug; see gh-16718.  When we can require Cython 3.0, y_local and
    # rank_local can be removed, and the closure weigh() can refer directly
    # to y and rank.
    cdef const ordered[:] y_local = y
    cdef intp_t i, first
    cdef float64_t t, u, v, w, s, sq
    cdef int64_t n = np.int64(len(x))
    cdef float64_t[::1] exchanges_weight = np.zeros(1, dtype=np.float64)
    # initial sort on values of x and, if tied, on values of y
    cdef intp_t[::1] perm = np.lexsort((y, x))
    cdef intp_t[::1] temp = np.empty(n, dtype=np.intp) # support structure

    if weigher is None:
        weigher = lambda x: 1./(1 + x)

    if rank is None:
        # To generate a rank array, we must first reverse the permutation
        # (to get higher ranks first) and then invert it.
        rank = np.empty(n, dtype=np.intp)
        rank[...] = perm[::-1]
        _invert_in_place(rank)

    cdef intp_t[:] rank_local = rank

    # weigh joint ties
    first = 0
    t = 0
    w = weigher(rank[perm[first]])
    s = w
    sq = w * w

    for i in range(1, n):
        if x[perm[first]] != x[perm[i]] or y[perm[first]] != y[perm[i]]:
            t += s * (i - first - 1) if additive else (s * s - sq) / 2
            first = i
            s = sq = 0

        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    t += s * (n - first - 1) if additive else (s * s - sq) / 2

    # weigh ties in x
    first = 0
    u = 0
    w = weigher(rank[perm[first]])
    s = w
    sq = w * w

    for i in range(1, n):
        if x[perm[first]] != x[perm[i]]:
            u += s * (i - first - 1) if additive else (s * s - sq) / 2
            first = i
            s = sq = 0

        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    u += s * (n - first - 1) if additive else (s * s - sq) / 2
    if first == 0: # x is constant (all ties)
        return np.nan

    # this closure recursively sorts sections of perm[] by comparing
    # elements of y[perm[]] using temp[] as support

    def weigh(intp_t offset, intp_t length):
        cdef intp_t length0, length1, middle, i, j, k
        cdef float64_t weight, residual

        if length == 1:
            return weigher(rank_local[perm[offset]])
        length0 = length // 2
        length1 = length - length0
        middle = offset + length0
        residual = weigh(offset, length0)
        weight = weigh(middle, length1) + residual
        if y_local[perm[middle - 1]] < y_local[perm[middle]]:
            return weight

        # merging
        i = j = k = 0

        while j < length0 and k < length1:
            if y_local[perm[offset + j]] <= y_local[perm[middle + k]]:
                temp[i] = perm[offset + j]
                residual -= weigher(rank_local[temp[i]])
                j += 1
            else:
                temp[i] = perm[middle + k]
                exchanges_weight[0] += weigher(rank_local[temp[i]]) * (
                    length0 - j) + residual if additive else weigher(
                    rank_local[temp[i]]) * residual
                k += 1
            i += 1

        perm[offset+i:offset+i+length0-j] = perm[offset+j:offset+length0]
        perm[offset:offset+i] = temp[0:i]
        return weight

    # weigh discordances
    weigh(0, n)

    # weigh ties in y
    first = 0
    v = 0
    w = weigher(rank[perm[first]])
    s = w
    sq = w * w

    for i in range(1, n):
        if y[perm[first]] != y[perm[i]]:
            v += s * (i - first - 1) if additive else (s * s - sq) / 2
            first = i
            s = sq = 0

        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    v += s * (n - first - 1) if additive else (s * s - sq) / 2
    if first == 0: # y is constant (all ties)
        return np.nan

    # weigh all pairs
    s = sq = 0
    for i in range(n):
        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    tot = s * (n - 1) if additive else (s * s - sq) / 2

    tau = ((tot - (v + u - t)) - 2. * exchanges_weight[0]
           ) / np.sqrt(tot - u) / np.sqrt(tot - v)
    return min(1., max(-1., tau))


# FROM MGCPY: https://github.com/neurodata/mgcpy

# Distance transforms used for MGC and Dcorr

# Columnwise ranking of data
@cython.wraparound(False)
@cython.boundscheck(False)
cdef _dense_rank_data(ndarray x):
    _, v = np.unique(x, return_inverse=True)
    return v + 1


@cython.wraparound(False)
@cython.boundscheck(False)
def _rank_distance_matrix(distx):
    # faster than np.apply_along_axis
    return np.hstack([_dense_rank_data(distx[:, i]).reshape(-1, 1) for i in range(distx.shape[0])])


@cython.wraparound(False)
@cython.boundscheck(False)
def _center_distance_matrix(distx, global_corr='mgc', is_ranked=True):
    cdef int n = distx.shape[0]
    cdef int m = distx.shape[1]
    cdef ndarray rank_distx = np.zeros(n * m)

    if is_ranked:
        rank_distx = _rank_distance_matrix(distx)

    if global_corr == "rank":
        distx = rank_distx.astype(np.float64, copy=False)

    # 'mgc' distance transform (col-wise mean) - default
    cdef ndarray exp_distx = np.repeat(((distx.mean(axis=0) * n) / (n-1)), n).reshape(-1, n).T

    # center the distance matrix
    cdef ndarray cent_distx = distx - exp_distx

    if global_corr != "mantel" and global_corr != "biased":
        np.fill_diagonal(cent_distx, 0)

    return cent_distx, rank_distx



# Centers each distance matrix and rank matrix
@cython.wraparound(False)
@cython.boundscheck(False)
def _transform_distance_matrix(distx, disty, global_corr='mgc', is_ranked=True):
    if global_corr == "rank":
        is_ranked = True

    cent_distx, rank_distx = _center_distance_matrix(distx, global_corr, is_ranked)
    cent_disty, rank_disty = _center_distance_matrix(disty, global_corr, is_ranked)

    transform_dist = {"cent_distx": cent_distx, "cent_disty": cent_disty,
                      "rank_distx": rank_distx, "rank_disty": rank_disty}

    return transform_dist


# MGC specific functions
@cython.wraparound(False)
@cython.boundscheck(False)
cdef _expected_covar(const float64_t[:, :] distx, const float64_t[:, :] disty,
                     const int64_t[:, :] rank_distx, const int64_t[:, :] rank_disty,
                     float64_t[:, :] cov_xy, float64_t[:] expectx,
                     float64_t[:] expecty):
    # summing up the element-wise product of A and B based on the ranks,
    # yields the local family of covariances
    cdef intp_t n = distx.shape[0]
    cdef float64_t a, b
    cdef intp_t i, j, k, l
    for i in range(n):
        for j in range(n):
            a = distx[i, j]
            b = disty[i, j]
            k = rank_distx[i, j]
            l = rank_disty[i, j]

            cov_xy[k, l] += a * b

            expectx[k] += a
            expecty[l] += b

    return np.asarray(expectx), np.asarray(expecty)


@cython.wraparound(False)
@cython.boundscheck(False)
cdef _covar_map(float64_t[:, :] cov_xy, intp_t nx, intp_t ny):
    # get covariances for each k and l
    cdef intp_t k, l
    for k in range(nx - 1):
        for l in range(ny - 1):
            cov_xy[k+1, l+1] += (cov_xy[k+1, l] + cov_xy[k, l+1] - cov_xy[k, l])

    return np.asarray(cov_xy)


@cython.wraparound(False)
@cython.boundscheck(False)
def _local_covariance(distx, disty, rank_distx, rank_disty):
    # convert float32 numpy array to int, as it will be used as array indices
    # [0 to n-1]
    rank_distx = np.asarray(rank_distx, np.int64) - 1
    rank_disty = np.asarray(rank_disty, np.int64) - 1

    cdef intp_t n = distx.shape[0]
    cdef intp_t nx = np.max(rank_distx) + 1
    cdef intp_t ny = np.max(rank_disty) + 1
    cdef ndarray cov_xy = np.zeros((nx, ny))
    cdef ndarray expectx = np.zeros(nx)
    cdef ndarray expecty = np.zeros(ny)

    # summing up the element-wise product of A and B based on the ranks,
    # yields the local family of covariances
    expectx, expecty = _expected_covar(distx, disty, rank_distx, rank_disty,
                                       cov_xy, expectx, expecty)

    cov_xy[:, 0] = np.cumsum(cov_xy[:, 0])
    expectx = np.cumsum(expectx)

    cov_xy[0, :] = np.cumsum(cov_xy[0, :])
    expecty = np.cumsum(expecty)

    cov_xy = _covar_map(cov_xy, nx, ny)
    # centering the covariances
    cov_xy = cov_xy - ((expectx.reshape(-1, 1) @ expecty.reshape(-1, 1).T) / n**2)

    return cov_xy


@cython.wraparound(False)
@cython.boundscheck(False)
def _local_correlations(distx, disty, global_corr='mgc'):
    transformed = _transform_distance_matrix(distx, disty, global_corr)

    # compute all local covariances
    cdef ndarray cov_mat = _local_covariance(
        transformed["cent_distx"],
        transformed["cent_disty"].T,
        transformed["rank_distx"],
        transformed["rank_disty"].T)

    # compute local variances for data A
    cdef ndarray local_varx = _local_covariance(
        transformed["cent_distx"],
        transformed["cent_distx"].T,
        transformed["rank_distx"],
        transformed["rank_distx"].T)
    local_varx = local_varx.diagonal()

    # compute local variances for data B
    cdef ndarray local_vary = _local_covariance(
        transformed["cent_disty"],
        transformed["cent_disty"].T,
        transformed["rank_disty"],
        transformed["rank_disty"].T)
    local_vary = local_vary.diagonal()

    # normalizing the covariances yields the local family of correlations

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        corr_mat = cov_mat / np.sqrt(local_varx.reshape(-1, 1) @ local_vary.reshape(-1, 1).T).real
        # avoid computational issues that may cause a few local correlations
        # to be negligibly larger than 1
        corr_mat[corr_mat > 1] = 1

    # set any local correlation to 0 if any corresponding local variance is
    # less than or equal to 0
    corr_mat[local_varx <= 0, :] = 0
    corr_mat[:, local_vary <= 0] = 0

    return corr_mat


cpdef double geninvgauss_logpdf(double x, double p, double b) noexcept nogil:
    return _geninvgauss_logpdf_kernel(x, p, b)


cdef double _geninvgauss_logpdf_kernel(double x, double p, double b) noexcept nogil:
    cdef double z, c

    if x <= 0:
        return -INFINITY

    z = cs.kve(p, b)
    if math.isinf(z):
        return NAN

    c = -math.log(2) - math.log(z) + b
    return c + (p - 1)*math.log(x) - b*(x + 1/x)/2


cdef double _geninvgauss_pdf(double x, void *user_data) noexcept nogil:
    # destined to be used in a LowLevelCallable
    cdef double p, b

    if x <= 0:
        return 0.

    p = (<double *>user_data)[0]
    b = (<double *>user_data)[1]

    return math.exp(_geninvgauss_logpdf_kernel(x, p, b))


cdef double _phi(double z) noexcept nogil:
    """evaluates the normal PDF. Used in `studentized_range`"""
    cdef double inv_sqrt_2pi = 0.3989422804014327
    return inv_sqrt_2pi * math.exp(-0.5 * z * z)


cdef double _logphi(double z) noexcept nogil:
    """evaluates the log of the normal PDF. Used in `studentized_range`"""
    cdef double log_inv_sqrt_2pi = -0.9189385332046727
    return log_inv_sqrt_2pi - 0.5 * z * z


cdef double _Phi(double z) noexcept nogil:
    """evaluates the normal CDF. Used in `studentized_range`"""
    # use a custom function because using cs.ndtr results in incorrect PDF at
    # q=0 on 32bit systems. Use a hardcoded 1/sqrt(2) constant rather than
    # math constants because they're not available on all systems.
    cdef double inv_sqrt_2 = 0.7071067811865475
    return 0.5 * math.erfc(-z * inv_sqrt_2)


cpdef double _studentized_range_cdf_logconst(double k, double df) noexcept:
    """Evaluates log of constant terms in the cdf integrand"""
    cdef double log_2 = 0.6931471805599453
    return (math.log(k) + (df / 2) * math.log(df)
            - (math.lgamma(df / 2) + (df / 2 - 1) * log_2))


cpdef double _studentized_range_pdf_logconst(double k, double df) noexcept:
    """Evaluates log of constant terms in the pdf integrand"""
    cdef double log_2 = 0.6931471805599453
    return (math.log(k) + math.log(k - 1) + (df / 2) * math.log(df)
            - (math.lgamma(df / 2) + (df / 2 - 1) * log_2))


cdef double _studentized_range_cdf(int n, double[2] integration_var,
                                   void *user_data) noexcept nogil:
    # evaluates the integrand of Equation (3) by Batista, et al [2]
    # destined to be used in a LowLevelCallable
    q = (<double *> user_data)[0]
    k = (<double *> user_data)[1]
    df = (<double *> user_data)[2]
    log_cdf_const = (<double *> user_data)[3]

    s = integration_var[1]
    z = integration_var[0]

    # suitable terms are evaluated within logarithms to avoid under/overflows
    log_terms = (log_cdf_const
                 + (df - 1) * math.log(s)
                 - (df * s * s / 2)
                 + _logphi(z))

    # multiply remaining term outside of log because it can be 0
    return math.exp(log_terms) * math.pow(_Phi(z + q * s) - _Phi(z), k - 1)


cdef double _studentized_range_cdf_asymptotic(double z, void *user_data) noexcept nogil:
    # evaluates the integrand of equation (2) by Lund, Lund, page 205. [4]
    # destined to be used in a LowLevelCallable
    q = (<double *> user_data)[0]
    k = (<double *> user_data)[1]

    return k * _phi(z) * math.pow(_Phi(z + q) - _Phi(z), k - 1)


cdef double _studentized_range_pdf(int n, double[2] integration_var,
                                   void *user_data) noexcept nogil:
    # evaluates the integrand of equation (4) by Batista, et al [2]
    # destined to be used in a LowLevelCallable
    q = (<double *> user_data)[0]
    k = (<double *> user_data)[1]
    df = (<double *> user_data)[2]
    log_pdf_const = (<double *> user_data)[3]

    z = integration_var[0]
    s = integration_var[1]

    # suitable terms are evaluated within logarithms to avoid under/overflows
    log_terms = (log_pdf_const
                 + df * math.log(s)
                 - df * s * s / 2
                 + _logphi(z)
                 + _logphi(s * q + z))

    # multiply remaining term outside of log because it can be 0
    return math.exp(log_terms) * math.pow(_Phi(s * q + z) - _Phi(z), k - 2)


cdef double _studentized_range_pdf_asymptotic(double z, void *user_data) noexcept nogil:
    # evaluates the integrand of equation (2) by Lund, Lund, page 205. [4]
    # destined to be used in a LowLevelCallable
    q = (<double *> user_data)[0]
    k = (<double *> user_data)[1]

    return k * (k - 1) * _phi(z) * _phi(z + q) * math.pow(_Phi(z + q) - _Phi(z), k - 2)


cdef double _studentized_range_moment(int n, double[3] integration_var,
                                      void *user_data) noexcept nogil:
    # destined to be used in a LowLevelCallable
    K = (<double *> user_data)[0]  # the Kth moment to calc.
    k = (<double *> user_data)[1]
    df = (<double *> user_data)[2]
    log_pdf_const = (<double *> user_data)[3]

    # Pull outermost integration variable out to pass as q to PDF
    q = integration_var[2]

    cdef double pdf_data[4]
    pdf_data[0] = q
    pdf_data[1] = k
    pdf_data[2] = df
    pdf_data[3] = log_pdf_const

    return (math.pow(q, K) *
            _studentized_range_pdf(4, integration_var, pdf_data))


cpdef double genhyperbolic_pdf(double x, double p, double a, double b) noexcept nogil:
    return math.exp(_genhyperbolic_logpdf_kernel(x, p, a, b))


cdef double _genhyperbolic_pdf(double x, void *user_data) noexcept nogil:
    # destined to be used in a LowLevelCallable
    cdef double p, a, b

    p = (<double *>user_data)[0]
    a = (<double *>user_data)[1]
    b = (<double *>user_data)[2]

    return math.exp(_genhyperbolic_logpdf_kernel(x, p, a, b))


cpdef double genhyperbolic_logpdf(
        double x, double p, double a, double b
        ) noexcept nogil:
    return _genhyperbolic_logpdf_kernel(x, p, a, b)


# logpdf is always negative, so use positive exception value
cdef double _genhyperbolic_logpdf(double x, void *user_data) noexcept nogil:
    # destined to be used in a LowLevelCallable
    cdef double p, a, b

    p = (<double *>user_data)[0]
    a = (<double *>user_data)[1]
    b = (<double *>user_data)[2]

    return _genhyperbolic_logpdf_kernel(x, p, a, b)


cdef double _genhyperbolic_logpdf_kernel(
        double x, double p, double a, double b
        ) noexcept nogil:
    cdef double t1, t2, t3, t4, t5

    t1 = _log_norming_constant(p, a, b)
    t2 = math.sqrt(1.0 + x*x)
    t3 = (p - 0.5) * math.log(t2)
    t4 = math.log(cs.kve(p - 0.5, a * t2)) - a * t2
    t5 = b * x

    return t1 + t3 + t4 + t5


cdef double _log_norming_constant(double p, double a, double b) noexcept nogil:
    cdef double t1, t2, t3, t4, t5, t6

    t1 = (a + b)*(a - b)
    t2 = p * 0.5 * math.log(t1)
    t3 = 0.5 * math.log(2 * PI)
    t4 = (p - 0.5) * math.log(a)
    t5 = math.sqrt(t1)
    t6 = math.log(cs.kve(p, t5)) - t5

    return t2 - t3 - t4 - t6


ctypedef fused real:
    float
    double
    long double


@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline int gaussian_kernel_estimate_inner(
    const real[:, :] points_,  const real[:, :] values_, const real[:, :] xi_,
    real[:, :] estimate, const real[:, :] cho_cov,
    int n, int m, int d, int p,
) noexcept nogil:
    cdef:
        int i, j, k
        real residual, arg, norm

    # Evaluate the normalisation
    norm = math.pow((2 * PI), (- d / 2.))
    for i in range(d):
        norm /= cho_cov[i, i]

    for i in range(n):
        for j in range(m):
            arg = 0
            for k in range(d):
                residual = (points_[i, k] - xi_[j, k])
                arg += residual * residual

            arg = math.exp(-arg / 2.) * norm
            for k in range(p):
                estimate[j, k] += values_[i, k] * arg

    return 0


@cython.wraparound(False)
@cython.boundscheck(False)
def gaussian_kernel_estimate(points, values, xi, cho_cov, dtype,
                             real _=0):
    """
    Evaluate a multivariate Gaussian kernel estimate.

    Parameters
    ----------
    points : array_like with shape (n, d)
        Data points to estimate from in d dimensions.
    values : real[:, :] with shape (n, p)
        Multivariate values associated with the data points.
    xi : array_like with shape (m, d)
        Coordinates to evaluate the estimate at in d dimensions.
    cho_cov : array_like with shape (d, d)
        (Lower) Cholesky factor of the covariance.

    Returns
    -------
    estimate : double[:, :] with shape (m, p)
        Multivariate Gaussian kernel estimate evaluated at the input coordinates.
    """
    cdef:
        real[:, :] points_, xi_, values_, estimate, cho_cov_
        int n, d, m, p

    n = points.shape[0]
    d = points.shape[1]
    m = xi.shape[0]
    p = values.shape[1]

    if xi.shape[1] != d:
        raise ValueError("points and xi must have same trailing dim")
    if cho_cov.shape[0] != d or cho_cov.shape[1] != d:
        raise ValueError("Covariance matrix must match data dims")

    # Rescale the data
    cho_cov_ = cho_cov.astype(dtype, copy=False)
    points_ = np.asarray(solve_triangular(cho_cov, points.T, lower=True).T,
                         dtype=dtype)
    xi_ = np.asarray(solve_triangular(cho_cov, xi.T, lower=True).T,
                     dtype=dtype)
    values_ = values.astype(dtype, copy=False)

    # Create the result array and evaluate the weighted sum
    estimate = np.zeros((m, p), dtype)

    with nogil:
        gaussian_kernel_estimate_inner(points_, values_, xi_,
                                       estimate, cho_cov_, n, m, d, p)

    return np.asarray(estimate)


@cython.wraparound(False)
@cython.boundscheck(False)
cdef real logsumexp(real a, real b):
    cdef:
        real c
    c = max(a, b)
    return c + math.log(math.exp(a-c) + math.exp(b-c))


@cython.wraparound(False)
@cython.boundscheck(False)
def gaussian_kernel_estimate_log(points, values, xi, cho_cov, dtype, real _=0):
    """
    def gaussian_kernel_estimate_log(points, real[:, :] values, xi, cho_cov)

    Evaluate the log of the estimated pdf on a provided set of points.

    Parameters
    ----------
    points : array_like with shape (n, d)
        Data points to estimate from in ``d`` dimensions.
    values : real[:, :] with shape (n, p)
        Multivariate values associated with the data points.
    xi : array_like with shape (m, d)
        Coordinates to evaluate the estimate at in ``d`` dimensions.
    cho_cov : array_like with shape (d, d)
        (Lower) Cholesky factor of the covariance.

    Returns
    -------
    estimate : double[:, :] with shape (m, p)
        The log of the multivariate Gaussian kernel estimate evaluated at the
        input coordinates.
    """
    cdef:
        real[:, :] points_, xi_, values_, log_values_, estimate
        int i, j, k
        int n, d, m, p
        real arg, residual, log_norm

    n = points.shape[0]
    d = points.shape[1]
    m = xi.shape[0]
    p = values.shape[1]

    if xi.shape[1] != d:
        raise ValueError("points and xi must have same trailing dim")
    if cho_cov.shape[0] != d or cho_cov.shape[1] != d:
        raise ValueError("Covariance matrix must match data dims")

    # Rescale the data
    points_ = np.asarray(solve_triangular(cho_cov, points.T, lower=True).T,
                         dtype=dtype)
    xi_ = np.asarray(solve_triangular(cho_cov, xi.T, lower=True).T,
                     dtype=dtype)
    values_ = values.astype(dtype, copy=False)

    log_values_ = np.empty((n, p), dtype)
    for i in range(n):
        for k in range(p):
            log_values_[i, k] = math.log(values_[i, k])

    # Evaluate the normalisation
    log_norm = (- d / 2) * math.log(2 * PI)
    for i in range(d):
        log_norm -= math.log(cho_cov[i, i])

    # Create the result array and evaluate the weighted sum
    estimate = np.full((m, p), fill_value=-np.inf, dtype=dtype)
    for i in range(n):
        for j in range(m):
            arg = 0
            for k in range(d):
                residual = (points_[i, k] - xi_[j, k])
                arg += residual * residual

            arg = -arg / 2 + log_norm
            for k in range(p):
                estimate[j, k] = logsumexp(estimate[j, k],
                                           arg + log_values_[i, k])

    return np.asarray(estimate)
