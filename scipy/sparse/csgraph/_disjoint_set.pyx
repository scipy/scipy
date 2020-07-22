"""
Disjoint set data structure
"""

import numpy as np
cimport numpy as np
cimport cython


cdef class DisjointSet:
    """Disjoint set data structure [1]_.

    Parameters
    ----------
    n : int
        The number of elements in the set.

    Notes
    -----
    .. versionadded: 1.6.0


    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Disjoint-set_data_structure


    Examples
    --------
    >>> from scipy.sparse.csgraph import DisjointSet

    Initialize a disjoint set with 4 elements:

    >>> dis = DisjointSet(4)

    Merge some subsets:

    >>> dis.union(0, 1)
    True
    >>> dis.union(2, 3)
    True
    >>> dis.union(3, 3)
    False

    Find a root element:

    >>> dis.find(1)
    0

    Attributes
    ----------
    n : int
        The number of elements in the set.
    nc : int
        The number of components/subsets.

    Methods
    -------
    union
    find
    """
    cdef:
        readonly np.npy_intp n
        readonly np.npy_intp nc
        readonly np.npy_intp[:] _parents
        readonly np.npy_intp[:] _sizes

    def __init__(DisjointSet self, np.intp_t n):
        self.n = n
        self.nc = n
        self._sizes = np.ones(n, dtype=np.intp)
        self._parents = np.arange(n, dtype=np.intp)

    def find(DisjointSet self, np.intp_t x):
        """Find the root element of `x`.

        Parameters
        ----------
        x : int
            Input element.

        Returns
        -------
        root : int
            Root element of `x`.
        """
        cdef np.npy_intp parent
        parents = self._parents
        parent = parents[x]
        while parent != parents[parent]:
            parent = parents[parent]
        parents[x] = parent
        return parent

    def union(DisjointSet self, np.intp_t a, np.intp_t b):
        """Merge the subsets of `a` and `b`. The smaller subset (the child) is
        merged into the the larger subset (the parent). If the subsets are of
        equal size, the parent is determined by subset root with the smallest
        index.

        Parameters
        ----------
        a : int
            Input element.
        b : int
            Input element.

        Returns
        -------
        merged : bool
            `True` if `a` and `b` were in disjoint sets, `False` otherwise.
        """
        a = self.find(a)
        b = self.find(b)
        if a == b:
            return False

        if b < a:
            a, b = b, a
        sizes = self._sizes
        parents = self._parents
        if sizes[a] < sizes[b]:
            parents[a] = b
            sizes[b] += sizes[a]
        else:
            parents[b] = a
            sizes[a] += sizes[b]

        self.nc -= 1
        return True
