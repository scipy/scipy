"""
Disjoint set data structure
"""
import collections


class DisjointSet:
    """ Disjoint set data structure for incremental connectivity queries.

    .. versionadded: 1.6.0

    Attributes
    ----------
    n_elements : int
        The number of elements in the set.
    n_subsets : int
        The number of subsets.

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

    >>> disjoint_set = DisjointSet([1, 2, 3, 'a', 'b'])

    Merge some subsets:

    >>> disjoint_set.merge(1, 2)
    True
    >>> disjoint_set.merge(3, 'a')
    True
    >>> disjoint_set.merge('a', 'b')
    True
    >>> disjoint_set.merge('b', 'b')
    False

    Find root elements:

    >>> disjoint_set[2]
    1
    >>> disjoint_set['b']
    3

    Test connectivity:

    >>> disjoint_set.connected(1, 2)
    True
    >>> disjoint_set.connected(1, 'b')
    False

    List elements in disjoint set:

    >>> list(disjoint_set)
    [1, 2, 3, 'a', 'b']
    """
    def __init__(self, elements=None):
        self.n_elements = 0
        self.n_subsets = 0
        self._sizes = {}
        self._parents = {}
        self._indices = collections.OrderedDict()
        if elements is not None:
            for x in elements:
                self.add(x)

    def __iter__(self):
        """Returns an iterator of the elements in the disjoint set.

        Elements are ordered by insertion order.
        """
        return iter(self._indices)

    def __contains__(self, x):
        return x in self._indices

    def __getitem__(self, x):
        """Find the root element of `x`.

        Parameters
        ----------
        x : hashable object
            Input element.

        Returns
        -------
        root : hashable object
            Root element of `x`.
        """
        if x not in self._indices:
            raise KeyError(x)

        # find by "path halving"
        parents = self._parents
        while self._indices[x] != self._indices[parents[x]]:
            parents[x] = parents[parents[x]]
            x = parents[x]
        return x

    def add(self, x):
        """Add element `x` to disjoint set"""
        self._sizes[x] = 1
        self._parents[x] = x
        self._indices[x] = len(self._indices)
        self.n_elements += 1
        self.n_subsets += 1

    def merge(self, x, y):
        """Merge the subsets of `x` and `y`.

        The smaller subset (the child) is merged into the larger subset (the
        parent). If the subsets are of equal size, the root element which was
        first inserted into the disjoint set is selected as the parent.

        Parameters
        ----------
        x, y : hashable object
            Elements to merge.

        Returns
        -------
        merged : bool
            True if `x` and `y` were in disjoint sets, False otherwise.
        """
        x = self[x]
        y = self[y]
        if self._indices[x] == self._indices[y]:
            return False

        sizes = self._sizes
        if (sizes[x], self._indices[y]) < (sizes[y], self._indices[x]):
            x, y = y, x
        self._parents[y] = x
        self._sizes[x] += self._sizes[y]
        self.n_subsets -= 1
        return True

    def connected(self, x, y):
        """Test whether `x` and `y` are in the same subset.

        Parameters
        ----------
        x, y : hashable object
            Elements to test.

        Returns
        -------
        result : bool
            True if `x` and `y` are in the same set, False otherwise.
        """
        x = self[x]
        y = self[y]
        return self._indices[x] == self._indices[y]
