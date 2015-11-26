"""
Spherical Voronoi Code

.. versionadded:: 0.17.0

"""
#
# Copyright (C)  Tyler Reddy, Ross Hemsley, Edd Edmondson,
#                Nikolai Nowaczyk, Joe Pitt-Francis, 2015.
#
# Distributed under the same BSD license as Scipy.
#

import numpy as np
import numpy.matlib
import scipy
import itertools
from scipy._lib._version import NumpyVersion

# Whether Numpy has stacked matrix linear algebra
HAS_NUMPY_VEC_DET = (NumpyVersion(np.__version__) >= '1.8.0')

__all__ = ['SphericalVoronoi']

def calc_circumcenters(tetrahedrons):
    """ Calculates the cirumcenters of the circumspheres of tetrahedrons.

    An implementation based on
    http://mathworld.wolfram.com/Circumsphere.html

    Parameters
    ----------
    tetrahedrons : an array of shape (N, 4, 3)
        consisting of N tetrahedrons defined by 4 points in 3D

    Returns
    ----------
    circumcenters : an array of shape (N, 3)
        consisting of the N circumcenters of the tetrahedrons in 3D

    """

    num = tetrahedrons.shape[0]
    a = np.concatenate((tetrahedrons, np.ones((num, 4, 1))), axis=2)

    sums = np.sum(tetrahedrons ** 2, axis=2)
    d = np.concatenate((sums[:, :, np.newaxis], a), axis=2)

    dx = np.delete(d, 1, axis=2)
    dy = np.delete(d, 2, axis=2)
    dz = np.delete(d, 3, axis=2)

    if HAS_NUMPY_VEC_DET:
        dx = np.linalg.det(dx)
        dy = -np.linalg.det(dy)
        dz = np.linalg.det(dz)
        a = np.linalg.det(a)
    else:
        dx = np.array([np.linalg.det(m) for m in dx])
        dy = -np.array([np.linalg.det(m) for m in dy])
        dz = np.array([np.linalg.det(m) for m in dz])
        a = np.array([np.linalg.det(m) for m in a])

    nominator = np.vstack((dx, dy, dz))
    denominator = 2*a
    return (nominator / denominator).T


def project_to_sphere(points, center, radius):
    """
    Projects the elements of points onto the sphere defined
    by center and radius.

    Parameters
    ----------
    points : array of floats of shape (npoints, ndim)
             consisting of the points in a space of dimension ndim
    center : array of floats of shape (ndim,)
            the center of the sphere to project on
    radius : float
            the radius of the sphere to project on

    returns: array of floats of shape (npoints, ndim)
            the points projected onto the sphere
    """

    lengths = scipy.spatial.distance.cdist(points, np.array([center]))
    return (points - center) / lengths * radius + center


class SphericalVoronoi:
    """ Voronoi diagrams on the surface of a sphere.

    .. versionadded:: 0.17.0

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, 3)
        Coordinates of points to construct a spherical
        Voronoi diagram from
    radius : float, optional
        Radius of the sphere (Default: 1)
    center : ndarray of floats, shape (3,)
        Center of sphere (Default: origin)

    Attributes
    ----------
    points : double array of shape (npoints, 3)
            the points in 3D to generate the Voronoi diagram from
    radius : double
            radius of the sphere
            Default: None (forces estimation, which is less precise)
    center : double array of shape (3,)
            center of the sphere
            Default: None (assumes sphere is centered at origin)
    vertices : double array of shape (nvertices, 3)
            Voronoi vertices corresponding to points
    regions : list of list of integers of shape (npoints, _ )
            the n-th entry is a list consisting of the indices
            of the vertices belonging to the n-th point in points

    Notes
    ----------
    The spherical Voronoi diagram algorithm proceeds as follows. The Convex
    Hull of the input points (generators) is calculated, and is equivalent to
    their Delaunay triangulation on the surface of the sphere [Caroli]_.
    A 3D Delaunay tetrahedralization is obtained by including the origin of
    the coordinate system as the fourth vertex of each simplex of the Convex
    Hull. The circumcenters of all tetrahedra in the system are calculated and
    projected to the surface of the sphere, producing the Voronoi vertices.
    The Delaunay tetrahedralization neighbour information is then used to
    order the Voronoi region vertices around each generator. The latter
    approach is substantially less sensitive to floating point issues than
    angle-based methods of Voronoi region vertex sorting.

    The surface area of spherical polygons is calculated by decomposing them
    into triangles and using L'Huilier's Theorem to calculate the spherical
    excess of each triangle [Weisstein]_. The sum of the spherical excesses is
    multiplied by the square of the sphere radius to obtain the surface area
    of the spherical polygon. For nearly-degenerate spherical polygons an area
    of approximately 0 is returned by default, rather than attempting the
    unstable calculation.

    Empirical assessment of spherical Voronoi algorithm performance suggests
    quadratic time complexity (loglinear is optimal, but algorithms are more
    challenging to implement). The reconstitution of the surface area of the
    sphere, measured as the sum of the surface areas of all Voronoi regions,
    is closest to 100 % for larger (>> 10) numbers of generators.

    References
    ----------

    .. [Caroli] Caroli et al. Robust and Efficient Delaunay triangulations of
                points on or close to a sphere. Research Report RR-7004, 2009.
    .. [Weisstein] "L'Huilier's Theorem." From MathWorld -- A Wolfram Web
                Resource. http://mathworld.wolfram.com/LHuiliersTheorem.html

    See Also
    --------
    Voronoi : Conventional Voronoi diagrams in N dimensions.

    Examples
    --------

    >>> from matplotlib import colors
    >>> from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    >>> import matplotlib.pyplot as plt
    >>> from scipy.spatial import SphericalVoronoi
    >>> from mpl_toolkits.mplot3d import proj3d
    >>> # set input data
    >>> points = np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0],
    ...                    [0, 1, 0], [0, -1, 0], [-1, 0, 0], ])
    >>> center = np.array([0, 0, 0])
    >>> radius = 1
    >>> # calculate spherical Voronoi diagram
    >>> sv = SphericalVoronoi(points, radius, center)
    >>> # sort vertices (optional, helpful for plotting)
    >>> sv.sort_vertices_of_regions()
    >>> # generate plot
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> # plot the unit sphere for reference (optional)
    >>> u = np.linspace(0, 2 * np.pi, 100)
    >>> v = np.linspace(0, np.pi, 100)
    >>> x = np.outer(np.cos(u), np.sin(v))
    >>> y = np.outer(np.sin(u), np.sin(v))
    >>> z = np.outer(np.ones(np.size(u)), np.cos(v))
    >>> ax.plot_surface(x, y, z, color='y', alpha=0.1)
    >>> # plot generator points
    >>> ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')
    >>> # plot Voronoi vertices
    >>> ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2],
    ...                    c='g')
    >>> # indicate Voronoi regions (as Euclidean polygons)
    >>> for region in sv.regions:
    ...    random_color = colors.rgb2hex(np.random.rand(3))
    ...    polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
    ...    polygon.set_color(random_color)
    ...    ax.add_collection3d(polygon)
    >>> plt.show()

    """

    def __init__(self, points, radius=None, center=None):
        """
        Initializes the object and starts the computation of the Voronoi
        diagram.

        points : The generator points of the Voronoi diagram assumed to be
         all on the sphere with radius supplied by the radius parameter and
         center supplied by the center parameter.
        radius : The radius of the sphere. Will default to 1 if not supplied.
        center : The center of the sphere. Will default to the origin if not
         supplied.
        """

        self.points = points
        if np.any(center):
            self.center = center
        else:
            self.center = np.zeros(3)
        if radius:
            self.radius = radius
        else:
            self.radius = 1
        self.vertices = None
        self.regions = None
        self._tri = None
        self._calc_vertices_regions()

    def _calc_vertices_regions(self):
        """
        Calculates the Voronoi vertices and regions of the generators stored
        in self.points. The vertices will be stored in self.vertices and the
        regions in self.regions.

        This algorithm was discussed at PyData London 2015 by
        Tyler Reddy, Ross Hemsley and Nikolai Nowaczyk
        """

        # perform 3D Delaunay triangulation on data set
        # (here ConvexHull can also be used, and is faster)
        self._tri = scipy.spatial.ConvexHull(self.points)

        # add the center to each of the simplices in tri to get the same
        # tetrahedrons we'd have gotten from Delaunay tetrahedralization
        tetrahedrons = self._tri.points[self._tri.simplices]
        tetrahedrons = np.insert(
            tetrahedrons,
            3,
            np.array([self.center]),
            axis=1
        )

        # produce circumcenters of tetrahedrons from 3D Delaunay
        circumcenters = calc_circumcenters(tetrahedrons)

        # project tetrahedron circumcenters to the surface of the sphere
        self.vertices = project_to_sphere(
            circumcenters,
            self.center,
            self.radius
        )

        # calculate regions from triangulation
        generator_indices = np.arange(self.points.shape[0])
        filter_tuple = np.where((np.expand_dims(self._tri.simplices,
                                -1) == generator_indices).any(axis=1))

        list_tuples_associations = zip(filter_tuple[1],
                                       filter_tuple[0])

        list_tuples_associations = sorted(list_tuples_associations,
                                          key=lambda t: t[0])

        # group by generator indices to produce
        # unsorted regions in nested list
        groups = []
        for k, g in itertools.groupby(list_tuples_associations,
                                      lambda t: t[0]):
            groups.append([element[1] for element in list(g)])

        self.regions = groups

    def sort_vertices_of_regions(self):
        """
         For each region in regions, it sorts the indices of the Voronoi
         vertices such that the resulting points are in a clockwise or
         counterclockwise order around the generator point.

         This is done as follows: Recall that the n-th region in regions
         surrounds the n-th generator in points and that the k-th
         Voronoi vertex in vertices is the projected circumcenter of the
         tetrahedron obtained by the k-th triangle in _tri.simplices (and the
         origin). For each region n, we choose the first triangle (=Voronoi
         vertex) in _tri.simplices and a vertex of that triangle not equal to
         the center n. These determine a unique neighbor of that triangle,
         which is then chosen as the second triangle. The second triangle
         will have a unique vertex not equal to the current vertex or the
         center. This determines a unique neighbor of the second triangle,
         which is then chosen as the third triangle and so forth. We proceed
         through all the triangles (=Voronoi vertices) belonging to the
         generator in points and obtain a sorted version of the vertices
         of its surrounding region.
        """

        for n in range(0, len(self.regions)):
            remaining = self.regions[n][:]
            sorted_vertices = []
            current_simplex = remaining[0]
            current_vertex = [k for k in self._tri.simplices[current_simplex]
                              if k != n][0]
            remaining.remove(current_simplex)
            sorted_vertices.append(current_simplex)
            while remaining:
                current_simplex = [
                    s for s in remaining
                    if current_vertex in self._tri.simplices[s]
                    ][0]
                current_vertex = [
                    s for s in self._tri.simplices[current_simplex]
                    if s != n and s != current_vertex
                    ][0]
                remaining.remove(current_simplex)
                sorted_vertices.append(current_simplex)
            self.regions[n] = sorted_vertices
