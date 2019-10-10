from __future__ import absolute_import

from cpython cimport bool
from libc cimport math
cimport cython
cimport numpy as np
from numpy.math cimport PI
from numpy cimport ndarray, int64_t, float64_t, intp_t

import numpy as np
import scipy.stats, scipy.special


cdef double von_mises_cdf_series(double k, double x, unsigned int p):
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
        R = 1. / (2 * n / k + R)
        V = R * (sn / n + V)

    with cython.cdivision(True):
        return 0.5 + x / (2 * PI) + V / PI


DEF SQRT2_PI = 0.79788456080286535588  # sqrt(2/pi)


cdef von_mises_cdf_normalapprox(k, x):
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
    for n in xrange(len(perm)-1, -1, -1):
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
    cdef intp_t i, j = 0, l = len(x)
    cdef intp_t[::1] perm = np.argsort(x, kind='quicksort')
    # The type of this array must be one of the supported types
    cdef int64_t[::1] result = np.ndarray(l, dtype=np.int64)

    # Find nans, if any, and assign them the lowest value
    for i in xrange(l - 1, -1, -1):
        if not np.isnan(x[perm[i]]):
            break
        result[perm[i]] = 0

    if i < l - 1:
        j = 1
        l = i + 1

    for i in xrange(l - 1):
        result[perm[i]] = j
        if x[perm[i]] != x[perm[i + 1]]:
            j += 1

    result[perm[i + 1]] = j
    return np.array(result, dtype=np.int64)


@cython.wraparound(False)
@cython.boundscheck(False)
def _weightedrankedtau(ordered[:] x, ordered[:] y, intp_t[:] rank, weigher, bool additive):
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

    # weigh joint ties
    first = 0
    t = 0
    w = weigher(rank[perm[first]])
    s = w
    sq = w * w

    for i in xrange(1, n):
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

    for i in xrange(1, n):
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
            return weigher(rank[perm[offset]])
        length0 = length // 2
        length1 = length - length0
        middle = offset + length0
        residual = weigh(offset, length0)
        weight = weigh(middle, length1) + residual
        if y[perm[middle - 1]] < y[perm[middle]]:
            return weight

        # merging
        i = j = k = 0

        while j < length0 and k < length1:
            if y[perm[offset + j]] <= y[perm[middle + k]]:
                temp[i] = perm[offset + j]
                residual -= weigher(rank[temp[i]])
                j += 1
            else:
                temp[i] = perm[middle + k]
                exchanges_weight[0] += weigher(rank[temp[i]]) * (
                    length0 - j) + residual if additive else weigher(
                    rank[temp[i]]) * residual
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

    for i in xrange(1, n):
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
    for i in xrange(n):
        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    tot = s * (n - 1) if additive else (s * s - sq) / 2

    tau = ((tot - (v + u - t)) - 2. * exchanges_weight[0]
           ) / np.sqrt(tot - u) / np.sqrt(tot - v)
    return min(1., max(-1., tau))
