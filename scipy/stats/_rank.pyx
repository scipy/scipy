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


DEF METHOD_AVERAGE = 0
DEF METHOD_MIN     = 1
DEF METHOD_MAX     = 2
DEF METHOD_DENSE   = 3
DEF METHOD_ORDINAL = 4

_tie_method_map = dict(average=METHOD_AVERAGE,
                       min=METHOD_MIN,
                       max=METHOD_MAX,
                       dense=METHOD_DENSE,
                       ordinal=METHOD_ORDINAL)

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
cdef np.ndarray[np.float64_t, ndim=1] _rankdata_fused(np.ndarray[array_data_type, ndim=1] b,
                                                      int tie_method):
    cdef unsigned int i, j, n, isize, dupcount, inext
    cdef np.ndarray[np.intp_t, ndim=1] order
    cdef np.ndarray[np.float64_t, ndim=1] ranks
    cdef double tie_rank
    cdef int total_tie_count = 0

    n = b.size
    ranks = _np.empty((n,))

    if tie_method == METHOD_ORDINAL:
        order = _np.argsort(b, kind="mergesort").astype(_np.intp)
    else:
        order = _np.argsort(b).astype(_np.intp)

    with nogil:
        if tie_method == METHOD_ORDINAL:
            for i in xrange(n):
                ranks[order[i]] = i + 1
        else:
            dupcount = 0
            for i in xrange(n):
                inext = i + 1
                if i == n - 1 or b[order[i]] != b[order[inext]]:
                    if tie_method == METHOD_AVERAGE:
                        tie_rank = inext - 0.5 * dupcount
                    elif tie_method == METHOD_MIN:
                        tie_rank = inext - dupcount
                    elif tie_method == METHOD_MAX:
                        tie_rank = inext
                    elif tie_method == METHOD_DENSE:
                        tie_rank = inext - dupcount - total_tie_count
                        total_tie_count += dupcount
                    for j in xrange(i - dupcount, inext):
                        ranks[order[j]] = tie_rank
                    dupcount = 0
                else:
                    dupcount += 1

    return ranks


@cython.boundscheck(False)
@cython.cdivision(True)
def rankdata(a, method='average'):
    """
    rankdata(a, method='average')

    Assign ranks to data, dealing with ties appropriately.

    Ranks begin at 1.  The `method` argument controls how ranks are assigned
    to equal values.  See [1]_ for further discussion of ranking methods.

    Parameters
    ----------
    a : array_like
        The array of values to be ranked.  The array is first flattened.
    method : str, optional
        The method used to assign ranks to tied elements.
        The options are 'average', 'min', 'max', 'dense' and 'ordinal'.

        'average':
            The average of the ranks that would have been assigned to
            all the tied values is assigned to each value.
        'min':
            The minimum of the ranks that would have been assigned to all
            the tied values is assigned to each value.  (This is also
            referred to as "competition" ranking.)
        'max':
            The maximum of the ranks that would have been assigned to all
            the tied values is assigned to each value.
        'dense':
            Like 'min', but the rank of the next highest element is assigned
            the rank immediately after those assigned to the tied elements.
        'ordinal':
            All values are given a distinct rank, corresponding to the order
            that the values occur in `a`.

        The default is 'average'.

    Returns
    -------
    ranks : ndarray
         An array of length equal to the size of `a`, containing rank
         scores.

    Notes
    -----
    All floating point types are converted to numpy.float64 before ranking.
    This may result in spurious ties if an input array of floats has a wider
    data type than numpy.float64 (e.g. numpy.float128).

    References
    ----------
    .. [1] "Ranking", http://en.wikipedia.org/wiki/Ranking

    Examples
    --------
    >>> from scipy.stats import rankdata
    >>> rankdata([0, 2, 3, 2])
    array([ 1. ,  2.5,  4. ,  2.5])
    >>> rankdata([0, 2, 3, 2], method='min')
    array([ 1.,  2.,  4.,  2.])
    >>> rankdata([0, 2, 3, 2], method='max')
    array([ 1.,  3.,  4.,  3.])
    >>> rankdata([0, 2, 3, 2], method='dense')
    array([ 1.,  2.,  3.,  2.])
    >>> rankdata([0, 2, 3, 2], method='ordinal')
    array([ 1.,  2.,  4.,  3.])
    """
    cdef np.ndarray[np.int64_t, ndim=1] b_int64
    cdef np.ndarray[np.uint64_t, ndim=1] b_uint64
    cdef np.ndarray[np.float64_t, ndim=1] b_float64
    cdef np.ndarray[np.float64_t, ndim=1] ranks

    # Convert the input to a 1-d numpy array.
    cdef np.ndarray b = _np.ravel(_np.asarray(a))

    cdef int tie_method

    if method not in _tie_method_map:
        raise ValueError("unknown method %r" % (method, ))
    tie_method = _tie_method_map[method]

    if b.size == 0:
        return _np.array([], dtype=_np.float64)

    if _np.issubdtype(b.dtype, _np.unsignedinteger):
        # Any unsigned type is converted to np.uint64.
        b_uint64 = _np.asarray(b, dtype=_np.uint64)
        ranks = _rankdata_fused(b_uint64, tie_method)
    elif _np.issubdtype(b.dtype, _np.integer):
        # Any integer type is converted to np.int64.
        b_int64 = _np.asarray(b, dtype=_np.int64)
        ranks = _rankdata_fused(b_int64, tie_method)
    else:
        # Anything else is converted to np.float64.
        b_float64 = _np.asarray(b, dtype=_np.float64)
        ranks = _rankdata_fused(b_float64, tie_method)

    return ranks


@cython.boundscheck(False)
@cython.cdivision(True)
def tiecorrect(rankvals):
    """
    tiecorrect(rankvals)

    Tie correction factor for ties in the Mann-Whitney U and
    Kruskal-Wallis H tests.

    Parameters
    ----------
    rankvals : array_like
        A 1-D sequence of ranks.  Typically this will be the array
        returned by `stats.rankdata`.

    Returns
    -------
    factor : float
        Correction factor for U or H.

    See Also
    --------
    rankdata : Assign ranks to the data
    mannwhitneyu : Mann-Whitney rank test
    kruskal : Kruskal-Wallis H test

    References
    ----------
    .. [1] Siegel, S. (1956) Nonparametric Statistics for the Behavioral
           Sciences.  New York: McGraw-Hill.

    Examples
    --------
    >>> from scipy.stats import tiecorrect, rankdata
    >>> tiecorrect([1, 2.5, 2.5, 4])
    0.9
    >>> ranks = rankdata([1, 3, 2, 4, 5, 7, 2, 8, 4])
    >>> ranks
    array([ 1. ,  4. ,  2.5,  5.5,  7. ,  8. ,  2.5,  9. ,  5.5])
    >>> tiecorrect(ranks)
    0.9833333333333333

    """
    cdef np.ndarray[np.float64_t, ndim=1] ranks
    cdef unsigned int i, inext, n, nties, T
    cdef np.ndarray[np.intp_t, ndim=1] order
    cdef np.ndarray[np.float64_t, ndim=1] sorted_
    cdef double factor

    ranks = _np.asarray(rankvals).astype(_np.float64)

    n = ranks.size
    if n < 2:
        return 1.0

    order = _np.argsort(ranks).astype(_np.intp)
    sorted_ = _np.empty((n,))

    with nogil:
        for i in xrange(n):
            sorted_[i] = ranks[order[i]]

        T = 0
        i = 0
        while i < n - 1:
            inext = i + 1
            if sorted_[i] == sorted_[inext]:
                nties = 1
                while i < n - 1 and sorted_[i] == sorted_[inext]:
                    nties += 1
                    i += 1
                    inext += 1
                T = T + nties**3 - nties
            i = inext

        factor = 1.0 - T / float(n**3 - n)

    return factor
