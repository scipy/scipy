"""
Disjoint set data structure
"""


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

    >>> dis.union(1, 2)
    True
    >>> dis.union(3, 'x')
    True
    >>> dis.union('x', 'y')
    True
    >>> dis.union('y', 'y')
    False

    Find root nodes:

    >>> dis[2]
    1
    >>> dis['y']
    3

    """
    def __init__(self):
        self.n_nodes = 0
        self.n_components = 0
        self._sizes = {}
        self._parents = {}
        self._indices = {}

    def __getitem__(self, x):
        """Find the root node of `x`.

        Parameters
        ----------
        x : hashable object
            Input node.

        Returns
        -------
        root : hashable object
            Root node of `x`.
        """
        if x not in self._parents:
            # add node
            self._sizes[x] = 1
            self._parents[x] = x
            self._indices[x] = len(self._indices)
            self.n_nodes += 1
            self.n_components += 1
            return x

        # find by "path halving"
        parents = self._parents
        while self._indices[x] != self._indices[parents[x]]:
            parents[x] = parents[parents[x]]
            x = parents[x]
        return x

    def union(self, a, b):
        """Merge the subsets of `a` and `b`.

        The smaller subset (the child) is merged into the larger subset (the
        parent).

        Parameters
        ----------
        a, b : hashable object
            Nodes to merge.

        Returns
        -------
        merged : bool
            `True` if `a` and `b` were in disjoint sets, `False` otherwise.
        """
        a = self[a]
        b = self[b]
        if self._indices[a] == self._indices[b]:
            return False

        if self._sizes[a] < self._sizes[b]:
            a, b = b, a
        self._parents[b] = a
        self._sizes[a] += self._sizes[b]
        self.n_components -= 1
        return True
