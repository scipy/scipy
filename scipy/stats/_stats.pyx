import warnings
import numpy as np
from numpy cimport int64_t, uint64_t
import scipy.special as special
from collections import namedtuple

cimport numpy as np
cimport cython

ctypedef fused ordered0:
    np.int32_t
    np.int64_t
    np.float32_t
    np.float64_t

ctypedef fused ordered1:
    np.int32_t
    np.int64_t
    np.float32_t
    np.float64_t

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline uint64_t pairs(uint64_t l):
    return (l * (l - 1)) // 2

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline double stat0(uint64_t l):
    return l * (l - 1.) * (2*l + 5)

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline double stat1(uint64_t l):
    return l * (l - 1.) * (l - 2)

@cython.wraparound(False)
@cython.boundscheck(False)
def _toranks(x):
    cdef long[:] perm = np.argsort(x, kind='quicksort')
    cdef long[:] rank = np.ndarray(len(perm), dtype=long)
    cdef uint64_t i, j = 0
    for i in xrange(len(x) - 1):
        rank[perm[i]] = j
        if x[perm[i]] != x[perm[i + 1]]:
            j += 1

    rank[perm[i + 1]] = j
    return rank

KendalltauResult = namedtuple('KendalltauResult', ('correlation', 'pvalue'))

@cython.wraparound(False)
@cython.boundscheck(False)
def _kendalltau(ordered0[:] x, ordered1[:] y):

    cdef uint64_t n = np.uint64(len(x))
    cdef long[:] temp = np.ndarray(n, dtype=long)  # support structure used by mergesort

    # this closure recursively sorts sections of perm[] by comparing
    # elements of y[perm[]] using temp[] as support
    # returns the number of swaps required by an equivalent bubble sort

    def mergesort(uint64_t offset, uint64_t length):
        cdef uint64_t exchcnt = 0, end, t, i, j, k, u
        # We use insertion sort on small arrays
        if length < 16:
            end = offset + length
            for i in xrange(offset + 1, end):
                t = perm[i]
                j = i
                u = perm[j - 1]
                while y[t] < y[u]:
                    exchcnt += 1
                    perm[j] = u
                    j -= 1
                    if offset == j:
                        break
                    u = perm[j - 1] 

                perm[j] = t
            return exchcnt

        cdef uint64_t length0 = length // 2
        cdef uint64_t length1 = length - length0
        cdef uint64_t middle = offset + length0
        exchcnt += mergesort(offset, length0)
        exchcnt += mergesort(middle, length1)
        if y[perm[middle - 1]] < y[perm[middle]]:
            return exchcnt

        # merging
        i = j = k = 0
        while j < length0 and k < length1:
            if y[perm[offset + j]] <= y[perm[middle + k]]:
                temp[i] = perm[offset + j]
                j += 1
            else:
                temp[i] = perm[middle + k]
                k += 1
                exchcnt += length0 - j
            i += 1

        perm[offset+i:offset+i+length0-j] = perm[offset+j:offset+length0]
        perm[offset:offset+i] = temp[0:i]
        return exchcnt

    # initial sort on values of x and, if tied, on values of y
    cdef long[:] perm = np.lexsort((y, x))
    
    # compute joint ties
    cdef uint64_t first = 0
    cdef uint64_t t = 0
    cdef uint64_t i
    for i in xrange(1, n):
        if x[perm[first]] != x[perm[i]] or y[perm[first]] != y[perm[i]]:
            t += pairs(i - first)
            first = i
    t += pairs(n - first)

    # compute ties in x
    first = 0
    cdef uint64_t u = 0
    cdef double u0 = 0, u1 = 0
    
    for i in xrange(1, n):
        if x[perm[first]] != x[perm[i]]:
            u += pairs(i - first)
            u0 += stat0(i - first)
            u1 += stat1(i - first)
            first = i

    if first == 0:  # All ties
        return KendalltauResult(np.nan, np.nan)
    u += pairs(n - first)
    u0 += stat0(n - first)
    u1 += stat1(n - first)

    # count exchanges
    cdef uint64_t exchanges = mergesort(0, n)

    # compute ties in y after mergesort with counting
    first = 0
    cdef uint64_t v = 0
    cdef double v0 = 0, v1 = 0
    for i in xrange(1, n):
        if y[perm[first]] != y[perm[i]]:
            v += pairs(i - first)
            v0 += stat0(i - first)
            v1 += stat1(i - first)
            first = i

    if first == 0:  # All ties
        return KendalltauResult(np.nan, np.nan)
    v += pairs(n - first)
    v0 += stat0(n - first)
    v1 += stat1(n - first)

    tot = (n * (n - 1)) // 2

    # Limit range to fix computational errors
    tau = min(1., max(-1., (tot - u - v + t - 2 * exchanges) / np.sqrt(tot - u) / np.sqrt(tot - v)))
    # (tot - u - v + t - 2 * exchanges) is approximately normally distributed with this variance
    var = (n * (n - 1) * (2*n + 5) - u0 - v0) / 18. + float(
        2 * u * v) / (n * (n - 1)) + float(u1 * v1) / (9 * n * (n - 1) * (n - 2))
    z = (tot - u - v + t - 2 * exchanges) / np.sqrt(var)
    prob = special.erfc(np.abs(z)/np.sqrt(2))
    return KendalltauResult(tau, prob)
