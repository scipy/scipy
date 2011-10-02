"""
Some functions related to ranking data, implemented in Cython for speed.
This file uses fused types, so Cython version 0.16 or later is required.

This module provides the following functions:
    rankdata
        Ranks the data, dealing with ties appropriately.
    tiecorrect
        Tie correction factor for the Mann-Whitney U test (stats.mannwhitnyu)
        and the Kruskal-Wallis H test (stats.kruskal).

"""

import numpy as _np
cimport numpy as np
cimport cython


ctypedef fused array_data_type:
    np.int64_t
    np.uint64_t
    np.float64_t


# The function _rankdata_fused uses the fused data type, array_data_type,
# so that arrays of unsigned integers, signed integers, floats are handled
# separately.  An alternative would be to convert any input array to, say,
# float64, but that would result in spurious ties when given large integers.
# For example float(2**60) == float(2**60 + 1), so the values 2**60 and
# 2**60 + 1 would be assigned the same rank if the input array was first
# converted to floats.

@cython.boundscheck(False)
@cython.cdivision(True)
cdef _rankdata_fused(np.ndarray[array_data_type, ndim=1] b):
    cdef unsigned int i, j, n, sumranks, dupcount, inext
    cdef np.ndarray[np.int_t, ndim=1] ivec
    cdef np.ndarray[array_data_type, ndim=1] svec
    cdef np.ndarray[np.float64_t, ndim=1] ranks
    cdef double averank

    n = b.size
    ivec = _np.argsort(b)
    svec = b[ivec]
    sumranks = 0
    dupcount = 0
    ranks = _np.zeros(n, float)
    for i in xrange(n):
        sumranks += i
        dupcount += 1
        inext = i + 1
        if i == n-1 or svec[i] != svec[inext]:
            averank = sumranks / float(dupcount) + 1
            for j in xrange(inext - dupcount, inext):
                ranks[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return ranks


@cython.boundscheck(False)
@cython.cdivision(True)
def rankdata(a):
    """
    rankata(a)

    Assign ranks to the data in `a`, dealing with ties appropriately.

    Equal values are assigned a rank that is the average of the ranks that
    would have been otherwise assigned to all of the values within that set.
    Ranks begin at 1.

    Parameters
    ----------
    a : array_like
        This array is first flattened.

    Returns
    -------
    ranks : ndarray
         An array of length equal to the size of `a`, containing rank scores.

    Notes
    -----
    All floating point types are converted to numpy.float64 before ranking.
    This may result in spurious ties if the input array has a wider data
    type than numpy.float64.

    Examples
    --------
    >>> rankdata([0, 2, 3, 2])
    array([ 1. ,  2.5,  4. ,  2.5])

    """
    cdef np.ndarray[np.int64_t, ndim=1] x_int64
    cdef np.ndarray[np.uint64_t, ndim=1] x_uint64
    cdef np.ndarray[np.float64_t, ndim=1] x_float64

    cdef np.ndarray b = _np.ravel(_np.asarray(a))

    if _np.issubdtype(b.dtype, _np.unsignedinteger):
        x_uint64 = b.astype(_np.uint64)
        ranks = _rankdata_fused(x_uint64)
    elif _np.issubdtype(b.dtype, _np.integer):
        x_int64 = b.astype(_np.int64)
        ranks = _rankdata_fused(x_int64)
    else:
        x_float64 = b.astype(_np.float64)
        ranks = _rankdata_fused(x_float64)

    return ranks


@cython.boundscheck(False)
@cython.cdivision(True)
def tiecorrect(np.ndarray[np.float64_t, ndim=1] rankvals):
    """
    tiecorrect(rankvals)

    Tie correction factor for ties in the Mann-Whitney U and
    Kruskal-Wallis H tests.

    Parameters
    ----------
    rankvals : 1-d ndarray of type np.float64
        For efficiency, this function requires a numpy array as its
        argument.  This is not a signficant inconvenience, because its
        primary use is in the mannwhitneyu and kruskal functions, where
        it is called with the result of stats.rankdata.

    Returns
    -------
    T : float
        Correction factor for U or H.

    See Also
    --------
    rankdata : Assign ranks to the data
    mannwhitney : Mann-Whitney rank test
    kruskal : Kruskal-Wallis H test

    Examples
    --------
    >>> ranks = rankdata([1, 3, 2, 4, 5, 7, 2, 8, 4])
    >>> ranks
    array([ 1. ,  4. ,  2.5,  5.5,  7. ,  8. ,  2.5,  9. ,  5.5])
    >>> tiecorrect(ranks)
    0.9833333333333333

    References
    ----------
    .. [1] Siegel, S. (1956) Nonparametric Statistics for the Behavioral
           Sciences.  New York: McGraw-Hill.

    """
    cdef unsigned int i, inext, n, nties, T
    cdef np.ndarray[np.int_t, ndim=1] posn
    cdef np.ndarray[np.float64_t, ndim=1] sorted
    cdef double t

    posn = _np.argsort(rankvals)
    sorted = rankvals[posn]
    n = len(sorted)
    if n < 2:
        return 1.0
    T = 0
    i = 0
    while i < n - 1:
        inext = i + 1
        if sorted[i] == sorted[inext]:
            nties = 1
            while i < n - 1 and sorted[i] == sorted[inext]:
                nties += 1
                i += 1
                inext += 1
            T = T + nties**3 - nties
        i = inext
    t = T / float(n**3 - n)
    return 1.0 - t
