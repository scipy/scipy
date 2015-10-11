"""
Some functions related to ranking data, implemented in Cython for speed.

This module provides the following functions:
    tiecorrect
        Tie correction factor for the Mann-Whitney U test (stats.mannwhitnyu)
        and the Kruskal-Wallis H test (stats.kruskal).

"""

import numpy as _np
cimport numpy as np
cimport cython


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
