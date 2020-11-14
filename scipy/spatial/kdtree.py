# Copyright Anne M. Archibald 2008
# Released under the scipy license
import numpy as np
import warnings
from .ckdtree import cKDTree, cKDTreeNode

__all__ = ['minkowski_distance_p', 'minkowski_distance',
           'distance_matrix',
           'Rectangle', 'KDTree']


def minkowski_distance_p(x, y, p=2):
    """
    Compute the pth power of the L**p distance between two arrays.

    For efficiency, this function computes the L**p distance but does
    not extract the pth root. If `p` is 1 or infinity, this is equal to
    the actual L**p distance.

    Parameters
    ----------
    x : (M, K) array_like
        Input array.
    y : (N, K) array_like
        Input array.
    p : float, 1 <= p <= infinity
        Which Minkowski p-norm to use.

    Examples
    --------
    >>> from scipy.spatial import minkowski_distance_p
    >>> minkowski_distance_p([[0,0],[0,0]], [[1,1],[0,1]])
    array([2, 1])

    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Find smallest common datatype with float64 (return type of this function) - addresses #10262.
    # Don't just cast to float64 for complex input case.
    common_datatype = np.promote_types(np.promote_types(x.dtype, y.dtype), 'float64')

    # Make sure x and y are NumPy arrays of correct datatype.
    x = x.astype(common_datatype)
    y = y.astype(common_datatype)

    if p == np.inf:
        return np.amax(np.abs(y-x), axis=-1)
    elif p == 1:
        return np.sum(np.abs(y-x), axis=-1)
    else:
        return np.sum(np.abs(y-x)**p, axis=-1)


def minkowski_distance(x, y, p=2):
    """
    Compute the L**p distance between two arrays.

    Parameters
    ----------
    x : (M, K) array_like
        Input array.
    y : (N, K) array_like
        Input array.
    p : float, 1 <= p <= infinity
        Which Minkowski p-norm to use.

    Examples
    --------
    >>> from scipy.spatial import minkowski_distance
    >>> minkowski_distance([[0,0],[0,0]], [[1,1],[0,1]])
    array([ 1.41421356,  1.        ])

    """
    x = np.asarray(x)
    y = np.asarray(y)
    if p == np.inf or p == 1:
        return minkowski_distance_p(x, y, p)
    else:
        return minkowski_distance_p(x, y, p)**(1./p)


class Rectangle(object):
    """Hyperrectangle class.

    Represents a Cartesian product of intervals.
    """
    def __init__(self, maxes, mins):
        """Construct a hyperrectangle."""
        self.maxes = np.maximum(maxes,mins).astype(float)
        self.mins = np.minimum(maxes,mins).astype(float)
        self.m, = self.maxes.shape

    def __repr__(self):
        return "<Rectangle %s>" % list(zip(self.mins, self.maxes))

    def volume(self):
        """Total volume."""
        return np.prod(self.maxes-self.mins)

    def split(self, d, split):
        """
        Produce two hyperrectangles by splitting.

        In general, if you need to compute maximum and minimum
        distances to the children, it can be done more efficiently
        by updating the maximum and minimum distances to the parent.

        Parameters
        ----------
        d : int
            Axis to split hyperrectangle along.
        split : float
            Position along axis `d` to split at.

        """
        mid = np.copy(self.maxes)
        mid[d] = split
        less = Rectangle(self.mins, mid)
        mid = np.copy(self.mins)
        mid[d] = split
        greater = Rectangle(mid, self.maxes)
        return less, greater

    def min_distance_point(self, x, p=2.):
        """
        Return the minimum distance between input and points in the hyperrectangle.

        Parameters
        ----------
        x : array_like
            Input.
        p : float, optional
            Input.

        """
        return minkowski_distance(0, np.maximum(0,np.maximum(self.mins-x,x-self.maxes)),p)

    def max_distance_point(self, x, p=2.):
        """
        Return the maximum distance between input and points in the hyperrectangle.

        Parameters
        ----------
        x : array_like
            Input array.
        p : float, optional
            Input.

        """
        return minkowski_distance(0, np.maximum(self.maxes-x,x-self.mins),p)

    def min_distance_rectangle(self, other, p=2.):
        """
        Compute the minimum distance between points in the two hyperrectangles.

        Parameters
        ----------
        other : hyperrectangle
            Input.
        p : float
            Input.

        """
        return minkowski_distance(0, np.maximum(0,np.maximum(self.mins-other.maxes,other.mins-self.maxes)),p)

    def max_distance_rectangle(self, other, p=2.):
        """
        Compute the maximum distance between points in the two hyperrectangles.

        Parameters
        ----------
        other : hyperrectangle
            Input.
        p : float, optional
            Input.

        """
        return minkowski_distance(0, np.maximum(self.maxes-other.mins,other.maxes-self.mins),p)


class KDTree(cKDTree):
    """
    kd-tree for quick nearest-neighbor lookup

    This class provides an index into a set of k-D points which
    can be used to rapidly look up the nearest neighbors of any point.

    Parameters
    ----------
    data : (N,K) array_like
        The data points to be indexed. This array is not copied, and
        so modifying this data will result in bogus results.
    leafsize : int, optional
        The number of points at which the algorithm switches over to
        brute-force.  Has to be positive.

    Notes
    -----
    The algorithm used is described in Maneewongvatana and Mount 1999.
    The general idea is that the kd-tree is a binary tree, each of whose
    nodes represents an axis-aligned hyperrectangle. Each node specifies
    an axis and splits the set of points based on whether their coordinate
    along that axis is greater than or less than a particular value.

    During construction, the axis and splitting point are chosen by the
    "sliding midpoint" rule, which ensures that the cells do not all
    become long and thin.

    The tree can be queried for the r closest neighbors of any given point
    (optionally returning only those within some maximum distance of the
    point). It can also be queried, with a substantial gain in efficiency,
    for the r approximate closest neighbors.

    For large dimensions (20 is already large) do not expect this to run
    significantly faster than brute force. High-dimensional nearest-neighbor
    queries are a substantial open problem in computer science.

    The tree also supports all-neighbors queries, both with arrays of points
    and with other kd-trees. These do use a reasonably efficient algorithm,
    but the kd-tree is not necessarily the best data structure for this
    sort of calculation.

    """

    class node:
        @staticmethod
        def _create(ckdtree_node=None):
            """Create either an inner or leaf node, wrapping a cKDTreeNode instance"""
            if ckdtree_node is None:
                return KDTree.node(ckdtree_node)
            elif ckdtree_node.split_dim == -1:
                return KDTree.leafnode(ckdtree_node)
            else:
                return KDTree.innernode(ckdtree_node)

        def __init__(self, ckdtree_node=None):
            if ckdtree_node is None:
                ckdtree_node = cKDTreeNode()
            self._node = ckdtree_node

        def __lt__(self, other):
            return id(self) < id(other)

        def __gt__(self, other):
            return id(self) > id(other)

        def __le__(self, other):
            return id(self) <= id(other)

        def __ge__(self, other):
            return id(self) >= id(other)

        def __eq__(self, other):
            return id(self) == id(other)

    class leafnode(node):
        @property
        def idx(self):
            return self._node.indices

        @property
        def children(self):
            return self._node.children

    class innernode(node):
        def __init__(self, ckdtreenode):
            assert isinstance(ckdtreenode, cKDTreeNode)
            super().__init__(ckdtreenode)
            self.less = KDTree.node._create(ckdtreenode.lesser)
            self.greater = KDTree.node._create(ckdtreenode.greater)

        @property
        def split_dim(self):
            return self._node.split_dim

        @property
        def split(self):
            return self._node.split

        @property
        def children(self):
            return self._node.children

    @property
    def tree(self):
        if not hasattr(self, "_tree"):
            self._tree = KDTree.node._create(super().tree)

        return self._tree

    def __init__(self, data, leafsize=10):
        data = np.asarray(data)
        if data.dtype.kind == 'c':
            raise TypeError("KDTree does not work with complex data")

        # Note KDTree has different default leafsize from cKDTree
        super().__init__(data, leafsize)

    def query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf):
        """
        Query the kd-tree for nearest neighbors

        Parameters
        ----------
        x : array_like, last dimension self.m
            An array of points to query.
        k : int, optional
            The number of nearest neighbors to return.
        eps : nonnegative float, optional
            Return approximate nearest neighbors; the kth returned value
            is guaranteed to be no further than (1+eps) times the
            distance to the real kth nearest neighbor.
        p : float, 1<=p<=infinity, optional
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        distance_upper_bound : nonnegative float, optional
            Return only neighbors within this distance. This is used to prune
            tree searches, so if you are doing a series of nearest-neighbor
            queries, it may help to supply the distance to the nearest neighbor
            of the most recent point.

        Returns
        -------
        d : float or array of floats
            The distances to the nearest neighbors.
            If x has shape tuple+(self.m,), then d has shape tuple if
            k is one, or tuple+(k,) if k is larger than one. Missing
            neighbors (e.g. when k > n or distance_upper_bound is
            given) are indicated with infinite distances.  If k is None,
            then d is an object array of shape tuple, containing lists
            of distances. In either case the hits are sorted by distance
            (nearest first).
        i : integer or array of integers
            The locations of the neighbors in self.data. i is the same
            shape as d.

        Examples
        --------
        >>> from scipy import spatial
        >>> x, y = np.mgrid[0:5, 2:8]
        >>> tree = spatial.KDTree(list(zip(x.ravel(), y.ravel())))
        >>> tree.data
        array([[0, 2],
               [0, 3],
               [0, 4],
               [0, 5],
               [0, 6],
               [0, 7],
               [1, 2],
               [1, 3],
               [1, 4],
               [1, 5],
               [1, 6],
               [1, 7],
               [2, 2],
               [2, 3],
               [2, 4],
               [2, 5],
               [2, 6],
               [2, 7],
               [3, 2],
               [3, 3],
               [3, 4],
               [3, 5],
               [3, 6],
               [3, 7],
               [4, 2],
               [4, 3],
               [4, 4],
               [4, 5],
               [4, 6],
               [4, 7]])
        >>> pts = np.array([[0, 0], [2.1, 2.9]])
        >>> tree.query(pts)
        (array([ 2.        ,  0.14142136]), array([ 0, 13]))
        >>> tree.query(pts[0])
        (2.0, 0)

        """
        x = np.asarray(x)
        if x.dtype.kind == 'c':
            raise TypeError("KDTree does not work with complex data")

        if k is None:
            # k=None, return all neighbors
            warnings.warn(
                "KDTree.query with k=None is deprecated and will be removed "
                "in SciPy 1.8.0. Use KDTree.query_ball_point instead.",
                DeprecationWarning)

            # Convert index query to a lists of distance and index,
            # sorted by distance
            def inds_to_hits(point, neighbors):
                dist = minkowski_distance(point, self.data[neighbors], p)
                hits = sorted([(d, i) for d, i in zip(dist, neighbors)])
                return [d for d, i in hits], [i for d, i in hits]

            x = np.asarray(x, dtype=np.float64)
            inds = super().query_ball_point(x, distance_upper_bound, p, eps)

            if isinstance(inds, list):
                return inds_to_hits(x, inds)

            dists = np.empty_like(inds)
            for idx in np.ndindex(inds.shape):
                dists[idx], inds[idx] = inds_to_hits(x[idx], inds[idx])

            return dists, inds

        d, i = super().query(x, k, eps, p, distance_upper_bound)
        if isinstance(i, int):
            i = np.intp(i)
        return d, i

    def query_ball_point(self, x, r, p=2., eps=0):
        """Find all points within distance r of point(s) x.

        Parameters
        ----------
        x : array_like, shape tuple + (self.m,)
            The point or points to search for neighbors of.
        r : positive float
            The radius of points to return.
        p : float, optional
            Which Minkowski p-norm to use.  Should be in the range [1, inf].
        eps : nonnegative float, optional
            Approximate search. Branches of the tree are not explored if their
            nearest points are further than ``r / (1 + eps)``, and branches are
            added in bulk if their furthest points are nearer than
            ``r * (1 + eps)``.

        Returns
        -------
        results : list or array of lists
            If `x` is a single point, returns a list of the indices of the
            neighbors of `x`. If `x` is an array of points, returns an object
            array of shape tuple containing lists of neighbors.

        Notes
        -----
        If you have many points whose neighbors you want to find, you may save
        substantial amounts of time by putting them in a KDTree and using
        query_ball_tree.

        Examples
        --------
        >>> from scipy import spatial
        >>> x, y = np.mgrid[0:5, 0:5]
        >>> points = np.c_[x.ravel(), y.ravel()]
        >>> tree = spatial.KDTree(points)
        >>> sorted(tree.query_ball_point([2, 0], 1))
        [5, 10, 11, 15]

        Query multiple points and plot the results:

        >>> import matplotlib.pyplot as plt
        >>> points = np.asarray(points)
        >>> plt.plot(points[:,0], points[:,1], '.')
        >>> for results in tree.query_ball_point(([2, 0], [3, 3]), 1):
        ...     nearby_points = points[results]
        ...     plt.plot(nearby_points[:,0], nearby_points[:,1], 'o')
        >>> plt.margins(0.1, 0.1)
        >>> plt.show()

        """
        x = np.asarray(x)
        if x.dtype.kind == 'c':
            raise TypeError("KDTree does not work with complex data")
        return super().query_ball_point(x, r, p, eps)

    def query_ball_tree(self, other, r, p=2., eps=0):
        """
        Find all pairs of points between `self` and `other` whose distance is at most r

        Parameters
        ----------
        other : KDTree instance
            The tree containing points to search against.
        r : float
            The maximum distance, has to be positive.
        p : float, optional
            Which Minkowski norm to use.  `p` has to meet the condition
            ``1 <= p <= infinity``.
        eps : float, optional
            Approximate search.  Branches of the tree are not explored
            if their nearest points are further than ``r/(1+eps)``, and
            branches are added in bulk if their furthest points are nearer
            than ``r * (1+eps)``.  `eps` has to be non-negative.

        Returns
        -------
        results : list of lists
            For each element ``self.data[i]`` of this tree, ``results[i]`` is a
            list of the indices of its neighbors in ``other.data``.

        Examples
        --------
        You can search all pairs of points between two kd-trees within a distance:

        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> from scipy.spatial import KDTree
        >>> np.random.seed(21701)
        >>> points1 = np.random.random((15, 2))
        >>> points2 = np.random.random((15, 2))
        >>> plt.figure(figsize=(6, 6))
        >>> plt.plot(points1[:, 0], points1[:, 1], "xk", markersize=14)
        >>> plt.plot(points2[:, 0], points2[:, 1], "og", markersize=14)
        >>> kd_tree1 = KDTree(points1)
        >>> kd_tree2 = KDTree(points2)
        >>> indexes = kd_tree1.query_ball_tree(kd_tree2, r=0.2)
        >>> for i in range(len(indexes)):
        ...     for j in indexes[i]:
        ...         plt.plot([points1[i, 0], points2[j, 0]],
        ...             [points1[i, 1], points2[j, 1]], "-r")
        >>> plt.show()

        """
        return super().query_ball_tree(other, r, p, eps)

    def query_pairs(self, r, p=2., eps=0):
        """
        Find all pairs of points in `self` whose distance is at most r.

        Parameters
        ----------
        r : positive float
            The maximum distance.
        p : float, optional
            Which Minkowski norm to use.  `p` has to meet the condition
            ``1 <= p <= infinity``.
        eps : float, optional
            Approximate search.  Branches of the tree are not explored
            if their nearest points are further than ``r/(1+eps)``, and
            branches are added in bulk if their furthest points are nearer
            than ``r * (1+eps)``.  `eps` has to be non-negative.

        Returns
        -------
        results : set
            Set of pairs ``(i,j)``, with ``i < j``, for which the corresponding
            positions are close.

        Examples
        --------
        You can search all pairs of points in a kd-tree within a distance:

        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> from scipy.spatial import KDTree
        >>> np.random.seed(21701)
        >>> points = np.random.random((20, 2))
        >>> plt.figure(figsize=(6, 6))
        >>> plt.plot(points[:, 0], points[:, 1], "xk", markersize=14)
        >>> kd_tree = KDTree(points)
        >>> pairs = kd_tree.query_pairs(r=0.2)
        >>> for (i, j) in pairs:
        ...     plt.plot([points[i, 0], points[j, 0]],
        ...             [points[i, 1], points[j, 1]], "-r")
        >>> plt.show()

        """
        return super().query_pairs(r, p, eps)

    def count_neighbors(self, other, r, p=2.):
        """
        Count how many nearby pairs can be formed.

        Count the number of pairs (x1,x2) can be formed, with x1 drawn
        from self and x2 drawn from ``other``, and where
        ``distance(x1, x2, p) <= r``.
        This is the "two-point correlation" described in Gray and Moore 2000,
        "N-body problems in statistical learning", and the code here is based
        on their algorithm.

        Parameters
        ----------
        other : KDTree instance
            The other tree to draw points from.
        r : float or one-dimensional array of floats
            The radius to produce a count for. Multiple radii are searched with
            a single tree traversal.
        p : float, 1<=p<=infinity, optional
            Which Minkowski p-norm to use

        Returns
        -------
        result : int or 1-D array of ints
            The number of pairs.

        Examples
        --------
        You can count neighbors number between two kd-trees within a distance:

        >>> import numpy as np
        >>> from scipy.spatial import KDTree
        >>> np.random.seed(21701)
        >>> points1 = np.random.random((5, 2))
        >>> points2 = np.random.random((5, 2))
        >>> kd_tree1 = KDTree(points1)
        >>> kd_tree2 = KDTree(points2)
        >>> kd_tree1.count_neighbors(kd_tree2, 0.2)
        9

        This number is same as the total pair number calculated by
        `query_ball_tree`:

        >>> indexes = kd_tree1.query_ball_tree(kd_tree2, r=0.2)
        >>> sum([len(i) for i in indexes])
        9

        """
        return super().count_neighbors(other, r, p)

    def sparse_distance_matrix(self, other, max_distance, p=2.):
        """
        Compute a sparse distance matrix

        Computes a distance matrix between two KDTrees, leaving as zero
        any distance greater than max_distance.

        Parameters
        ----------
        other : KDTree

        max_distance : positive float

        p : float, optional

        Returns
        -------
        result : dok_matrix
            Sparse matrix representing the results in "dictionary of keys" format.

        Examples
        --------
        You can compute a sparse distance matrix between two kd-trees:

        >>> import numpy as np
        >>> from scipy.spatial import KDTree
        >>> np.random.seed(21701)
        >>> points1 = np.random.random((5, 2))
        >>> points2 = np.random.random((5, 2))
        >>> kd_tree1 = KDTree(points1)
        >>> kd_tree2 = KDTree(points2)
        >>> sdm = kd_tree1.sparse_distance_matrix(kd_tree2, 0.3)
        >>> sdm.toarray()
        array([[0.20220215, 0.14538496, 0.,         0.10257199, 0.        ],
            [0.13491385, 0.27251306, 0.,         0.18793787, 0.        ],
            [0.19262396, 0.,         0.,         0.25795122, 0.        ],
            [0.14859639, 0.07076002, 0.,         0.04065851, 0.        ],
            [0.17308768, 0.,         0.,         0.24823138, 0.        ]])

        You can check distances above the `max_distance` are zeros:

        >>> from scipy.spatial import distance_matrix
        >>> distance_matrix(points1, points2)
        array([[0.20220215, 0.14538496, 0.43588092, 0.10257199, 0.4555495 ],
            [0.13491385, 0.27251306, 0.65944131, 0.18793787, 0.68184154],
            [0.19262396, 0.34121593, 0.72176889, 0.25795122, 0.74538858],
            [0.14859639, 0.07076002, 0.48505773, 0.04065851, 0.50043591],
            [0.17308768, 0.32837991, 0.72760803, 0.24823138, 0.75017239]])

        """
        return super().sparse_distance_matrix(other, max_distance, p)


def distance_matrix(x, y, p=2, threshold=1000000):
    """
    Compute the distance matrix.

    Returns the matrix of all pair-wise distances.

    Parameters
    ----------
    x : (M, K) array_like
        Matrix of M vectors in K dimensions.
    y : (N, K) array_like
        Matrix of N vectors in K dimensions.
    p : float, 1 <= p <= infinity
        Which Minkowski p-norm to use.
    threshold : positive int
        If ``M * N * K`` > `threshold`, algorithm uses a Python loop instead
        of large temporary arrays.

    Returns
    -------
    result : (M, N) ndarray
        Matrix containing the distance from every vector in `x` to every vector
        in `y`.

    Examples
    --------
    >>> from scipy.spatial import distance_matrix
    >>> distance_matrix([[0,0],[0,1]], [[1,0],[1,1]])
    array([[ 1.        ,  1.41421356],
           [ 1.41421356,  1.        ]])

    """

    x = np.asarray(x)
    m, k = x.shape
    y = np.asarray(y)
    n, kk = y.shape

    if k != kk:
        raise ValueError("x contains %d-dimensional vectors but y contains %d-dimensional vectors" % (k, kk))

    if m*n*k <= threshold:
        return minkowski_distance(x[:,np.newaxis,:],y[np.newaxis,:,:],p)
    else:
        result = np.empty((m,n),dtype=float)  # FIXME: figure out the best dtype
        if m < n:
            for i in range(m):
                result[i,:] = minkowski_distance(x[i],y,p)
        else:
            for j in range(n):
                result[:,j] = minkowski_distance(x,y[j],p)
        return result
