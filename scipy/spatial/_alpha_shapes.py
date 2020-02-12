"""
Alpha Shapes Code

.. versionadded:: 1.5.0

"""
import itertools
import numpy as np
from scipy.spatial import Delaunay
from . import _disjoint_set


__all__ = ['AlphaShapes']


class AlphaShapes:
    """
    AlphaShapes(points, qhull_options=None, rcond=1e-15)

    Alpha shapes in N dimensions.

    .. versionadded:: 1.5.0

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points to triangulate
    qhull_options : str, optional
        Additional options to pass to Qhull. See Qhull manual for
        details. Option "Qt" is always enabled.
        Default:"Qbb Qc Qz Qx Q12" for ndim > 4 and "Qbb Qc Qz Q12" otherwise.
        Incremental mode omits "Qz".
    rcond : float, optional
        Cutoff for small singular values for handling degenerate simplices.

    Attributes
    ----------
    points : ndarray of double, shape (npoints, ndim)
        Coordinates of input points.
    simplices : ndarray of ints, shape (nsimplex, ndim+1)
        Indices of the points forming the simplices in the alpha complex.
        For 2-D, the points are oriented counterclockwise.
    radii : ndarray of double, shape (nsimplex,)
        Circumradii of the simplices in sorted order.
    circumcenters : ndarray of double, shape (nsimplex, ndim)
        Circumcenters of the simplices.

    Raises
    ------
    QhullError
        Raised when Qhull encounters an error condition, such as
        geometrical degeneracy when options to resolve are not enabled.
    ValueError
        Raised if an incompatible array is given as input.

    Notes
    -----
    The tessellation is computed using the Qhull library
    `Qhull library <http://www.qhull.org/>`__.

    .. note::

       Unless you pass in the Qhull option "QJ", Qhull does not
       guarantee that each input point appears as a vertex in the
       Delaunay triangulation. Omitted points are listed in the
       `coplanar` attribute.

    Examples
    --------
    Alpha complex of a set of points:

    >>> import numpy as np
    >>> from scipy.spatial import AlphaShapes
    >>> rng = np.random.RandomState(seed=0)
    >>> points = rng.normal(size=(12, 2))
    >>> alpha = AlphaShapes(points)

    The simplices are ordered by circumradius:

    >>> alpha.radii
    array([ 0.41419783 0.5038057  0.51805934 0.52050488 0.53245625 0.55120833
            0.58252112 0.60508938 0.66602254 0.87210292 0.90357423 0.98832633
            1.02969871 1.41936012 1.770874   1.86824635 4.2987339 ])

    The number of connected components in the complex decreases as the
    threshold increases:

    >>> alpha.thresholds
    array([ 0.         1.41936012 1.02969871 0.60508938 0.58252112 0.58252112
            0.53245625 0.52050488 0.51805934 0.5038057  0.41419783 0.41419783])

    We can plot the alpha complex at a radius threshold which is large
    enough to include all points:

    >>> import matplotlib.pyplot as plt
    >>> from matplotlib import collections
    >>> facets = alpha.get_surface_facets(alpha.thresholds[1])
    >>> fig, ax = plt.subplots()
    >>> ax.scatter(points.T[0], points.T[1])
    >>> lc = collections.LineCollection(points[facets])
    >>> ax.add_collection(lc)
    >>> plt.show()
    """
    def __init__(self, points, qhull_options=None, rcond=1e-15):
        delaunay = Delaunay(points, qhull_options=qhull_options)
        self.points = delaunay.points
        self.simplices = delaunay.simplices

        self._calculate_circumcenters(rcond)
        self._calculate_circumradii()
        self._calculate_surface_intervals()
        self._calculate_connectivity_thresholds()

    def _calculate_circumcenters(self, rcond):
        tetrahedra = self.points[self.simplices]

        # build overdetermined system of linear equations (Ax=b)
        gramian = np.einsum('lij,lkj->lik', tetrahedra, tetrahedra)
        A = 2 * (gramian - np.roll(gramian, 1, axis=1))

        squared_norms = np.einsum('ij,ij->i', self.points, self.points)
        squared_tet_norms = squared_norms[self.simplices]
        b = squared_tet_norms - np.roll(squared_tet_norms, 1, axis=1)

        # handle rank deficiencies with Moore-Penrose pseudoinverse
        penrose = np.linalg.pinv(A, rcond=rcond)
        tp = np.einsum('lji,ljk->lik', tetrahedra, penrose)
        self.circumcenters = np.einsum('lij,lj->li', tp, b)

    def _calculate_circumradii(self):
        # calculate circumradii of each tetrahedron
        deltas = self.circumcenters - self.points[self.simplices[:, 0]]
        self.radii = np.linalg.norm(deltas, axis=1)

        # now sort the radii, simplices, and circumcenters by radius
        indices = np.argsort(self.radii)
        self.radii = self.radii[indices]
        self.circumcenters = self.circumcenters[indices]
        self.simplices = self.simplices[indices]

    def _calculate_surface_intervals(self):
        dim = self.points.shape[1]
        nsimplices = len(self.simplices)

        # decompose simplices into facets
        facets = []
        for indices in itertools.combinations(range(dim + 1), dim):
            facets.append(self.simplices[:, indices])
        facets = np.hstack(facets).reshape((nsimplices * (dim + 1), dim))
        facets = np.sort(facets)

        # calculate intervals in which facets are surface facets
        unique_facets, start, inverse = np.unique(facets, axis=0,
                                                  return_index=True,
                                                  return_inverse=True)
        _, end = np.unique(inverse[::-1], return_index=True)
        end = len(inverse) - 1 - end

        # facets which only appear once in list can only be surface facets
        indices = np.where(start == end)[0]
        end[indices] = nsimplices * (dim + 1)

        self._start = start
        self._end = end
        self._unique_facets = unique_facets

    def _calculate_connectivity_thresholds(self):
        n = len(self.points)
        uf = _disjoint_set.DisjointSet(n)
        self.thresholds = np.zeros(n)
        for ii, s in enumerate(self.simplices):
            for i, j in itertools.combinations(s, 2):
                if uf.merge(i, j):
                    n -= 1
                    self.thresholds[n] = self.radii[ii]

    def get_surface_facets(self, alpha):
        dim = self.points.shape[1]
        indices = np.where(self.radii <= alpha)[0]
        if not len(indices):
            return np.array([])

        index = (np.max(indices) + 1) * (dim + 1)
        indices = np.where((self._start < index) & (self._end >= index))[0]
        return self._unique_facets[indices]
