import warnings
import numpy as np
from numpy cimport int64_t, uint64_t
import scipy.special as special
from collections import namedtuple

cimport numpy as np
cimport cython

ctypedef fused ordered:
    char
    short
    int
    long
    float
    double

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline uint64_t pairs(uint64_t l):
    return (l * (l - 1)) // 2

KendalltauResult = namedtuple('KendalltauResult', ('correlation', 'pvalue'))

@cython.wraparound(False)
@cython.boundscheck(False)
def _kendalltau(ordered[:] x, ordered[:] y):

    cdef uint64_t n = np.uint64(len(x))
    cdef long[:] temp = np.ndarray(n, dtype=long) # support structure used by mergesort

    # this closure recursively sorts sections of perm[] by comparing
    # elements of y[perm[]] using temp[] as support
    # returns the number of swaps required by an equivalent bubble sort

    def mergesort(uint64_t offset, uint64_t length):
        cdef uint64_t exchcnt = 0, end, t, i, j, k, u
        # We use insertion sort on small arrays
        if length < 32:
            end = offset + length;
            for i in xrange(offset + 1, end):
                t = perm[i]
                j = i;
                u = perm[j - 1]
                while y[t] < y[u]:
                    exchcnt += 1
                    perm[j] = u
                    j -= 1
                    if offset == j:
                        break
                    u = perm[j - 1] 

                perm[j] = t
            return exchcnt;

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
                exchcnt += length0 - j;
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
    for i in xrange(1, n):
        if x[perm[first]] != x[perm[i]]:
            u += pairs(i - first)
            first = i
    u += pairs(n - first)

    # count exchanges
    cdef uint64_t exchanges = mergesort(0, n)
    # compute ties in y after mergesort with counting
    first = 0
    cdef uint64_t v = 0
    
    for i in xrange(1, n):
        if y[perm[first]] != y[perm[i]]:
            v += pairs(i - first)
            first = i
    v += pairs(n - first)

    tot = (n * (n - 1)) // 2
    if tot == u or tot == v:
        # Special case for all ties in one of the ranks
        return KendalltauResult(np.nan, np.nan)

    # Limit range to fix computational errors
    tau = min(1, max(-1, ((tot + t) - (v + u) - 2 * exchanges) / np.sqrt(tot - u) / np.sqrt(tot - v)))

    # what follows reproduces the ending of Gary Strangman's original
    # stats.kendalltau() in SciPy
    if u == 0 and v == 0:
        svar = (4.0 * n + 10.0) / (9.0 * n * (n - 1))
        z = tau / np.sqrt(svar)
        prob = special.erfc(np.abs(z) / 1.4142136)
        return KendalltauResult(tau, prob)
    
    return KendalltauResult(tau, np.nan)
