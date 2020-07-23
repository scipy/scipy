"""
Disjoint set data structure
"""

import numpy as np
cimport numpy as np
cimport cython


cdef class DisjointSet:
    """ Disjoint set data structure for incremental connectivity queries.

    .. versionadded: 1.6.0

    Parameters
    ----------
    n : int
        The number of elements in the set.

    Attributes
    ----------
    n : int
        The number of elements in the set.
    n_components : int
        The number of components/subsets.

    Methods
    -------
    union
    find

    Notes
    -----
    This class implements the disjoint set [1]_, also known as the *union-find*
    data structure. The *find* method implements the *path compression*
    variant. The *union* method implements the *union by size* variant.

    Element indices are integers in the range `[0, 1, ..., n - 1]`.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Disjoint-set_data_structure

    See Also
    --------
    connected_components : Analyze the connected components of a sparse graph.

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

    """
    cdef:
        readonly np.npy_intp n
        readonly np.npy_intp n_components
        readonly np.npy_intp[:] _parents
        readonly np.npy_intp[:] _sizes

    def __init__(DisjointSet self, np.intp_t n):
        self.n = n
        self.n_components = n
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
        """Merge the subsets of `a` and `b`.

        The smaller subset (the child) is merged into the the larger subset
        (the parent). If the subsets are of equal size, the parent is
        determined by subset root with the smallest index.

        Parameters
        ----------
        a, b : int
            Element indices to merge.

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

        self.n_components -= 1
        return True
