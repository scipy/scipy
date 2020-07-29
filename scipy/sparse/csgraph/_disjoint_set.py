"""
Disjoint set data structure
"""
import collections


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
    merge
    __getitem__

    Notes
    -----
    This class implements the disjoint set [1]_, also known as the *union-find*
    or *merge-find* data structure. The *find* operation (implemented in
    `__getitem__`)* implements the *path halving* variant. The *merge* method
    implements the *merge by size* variant.

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

    >>> disjoint_set = DisjointSet()

    Merge some subsets:

    >>> disjoint_set.merge(1, 2)
    True
    >>> disjoint_set.merge(3, 'a')
    True
    >>> disjoint_set.merge('a', 'b')
    True
    >>> disjoint_set.merge('b', 'b')
    False

    Find root nodes:

    >>> disjoint_set[2]
    1
    >>> disjoint_set['b']
    3

    """
    def __init__(self):
        self.n_nodes = 0
        self.n_components = 0
        self._sizes = {}
        self._parents = {}
        self._indices = collections.OrderedDict()

    def __iter__(self):
        """Returns an iterator of the nodes in the disjoint set.

        Nodes in the disjoint set consist of all inputs to the `__getitem__`,
        `merge`, and `connected` methods. Elements are ordered by insertion
        order.
        """
        return iter(self._indices)

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

    def merge(self, x, y):
        """Merge the subsets of `x` and `y`.

        The smaller subset (the child) is merged into the larger subset (the
        parent). If neither `x` nor `y` has been previously seen, `x` is
        inserted into the disjoint set before `y`.

        Parameters
        ----------
        x, y : hashable object
            Nodes to merge.

        Returns
        -------
        merged : bool
            True if `x` and `y` were in disjoint sets, False otherwise.
        """
        x = self[x]
        y = self[y]
        if self._indices[x] == self._indices[y]:
            return False

        if self._sizes[x] < self._sizes[y]:
            x, y = y, x
        self._parents[y] = x
        self._sizes[x] += self._sizes[y]
        self.n_components -= 1
        return True

    def connected(self, x, y):
        """Test whether `x` and `y` are in the same component/set.

        If neither `x` nor `y` has been previously seen, `x` is
        inserted into the disjoint set before `y`.

        Parameters
        ----------
        x, y : hashable object
            Nodes to test.

        Returns
        -------
        result : bool
            True if `x` and `y` are in the same set, False otherwise.
        """
        x = self[x]
        y = self[y]
        return self._indices[x] == self._indices[y]
