"""
Disjoint set data structure
"""


import numpy as np


class DisjointSet:
    """ Disjoint set data structure for incremental connectivity queries.

    .. versionadded: 1.6.0

    Attributes
    ----------
    n_nodes : int
        The number of nodes in the set.
    n_components : int
        The number of components/subsets.

    Methods
    -------
    union
    find

    Notes
    -----
    This class implements the disjoint set [1]_, also known as the *union-find*
    data structure. The *find* method implements the *path halving*
    variant. The *union* method implements the *union by size* variant.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Disjoint-set_data_structure

    See Also
    --------
    connected_components : Analyze the connected components of a sparse graph.

    Examples
    --------
    >>> from scipy.sparse.csgraph import DisjointSet

    Initialize a disjoint set:

    >>> dis = DisjointSet()

    Merge some subsets:

    >>> dis.union(0, 1)
    True
    >>> dis.union(2, 3)
    True
    >>> dis.union(3, 3)
    False

    Find a root node:

    >>> dis.find(1)
    0

    """
    def __init__(self):
        self.n_nodes = 0
        self.n_components = 0
        self._sizes = {}
        self._parents = {}

    def find(self, x):
        """Find the root node of `x`.

        Parameters
        ----------
        x : int
            Input node.

        Returns
        -------
        root : int
            Root node of `x`.
        """
        if not np.issubdtype(type(x), np.integer):
            raise TypeError("`x` must be an integer")

        if x not in self._parents:
            self._sizes[x] = 1
            self._parents[x] = x
            self.n_nodes += 1
            self.n_components += 1
            return x

        # find by "path halving"
        parents = self._parents
        while x != parents[x]:
            parents[x] = parents[parents[x]]
            x = parents[x]
        return x

    def union(self, a, b):
        """Merge the subsets of `a` and `b`.

        The smaller subset (the child) is merged into the larger subset (the
        parent). If the subsets are of equal size, the parent is determined by
        subset root with the smallest index.

        Parameters
        ----------
        a, b : int
            Node indices to merge.

        Returns
        -------
        merged : bool
            `True` if `a` and `b` were in disjoint sets, `False` otherwise.
        """
        if not np.issubdtype(type(a), np.integer):
            raise TypeError("`a` must be an integer")
        if not np.issubdtype(type(b), np.integer):
            raise TypeError("`b` must be an integer")

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
