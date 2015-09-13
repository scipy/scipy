.. _qhulltutorial:

Spatial data structures and algorithms (`scipy.spatial`)
========================================================

.. currentmodule:: scipy.spatial

`scipy.spatial` can compute triangulations, Voronoi diagrams, and
convex hulls of a set of points, by leveraging the `Qhull
<http://qhull.org/>`__ library.

Moreover, it contains `KDTree` implementations for nearest-neighbor point
queries, and utilities for distance computations in various metrics.

Delaunay triangulations
-----------------------

The Delaunay triangulation is a subdivision of a set of points into a
non-overlapping set of triangles, such that no point is inside the
circumcircle of any triangle. In practice, such triangulations tend to
avoid triangles with small angles.

Delaunay triangulation can be computed using `scipy.spatial` as follows:

.. plot::

   >>> from scipy.spatial import Delaunay
   >>> points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
   >>> tri = Delaunay(points)
   
   We can visualize it:
   
   >>> import matplotlib.pyplot as plt
   >>> plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
   >>> plt.plot(points[:,0], points[:,1], 'o')
   
   And add some further decorations:
   
   >>> for j, p in enumerate(points):
   ...     plt.text(p[0]-0.03, p[1]+0.03, j, ha='right') # label the points
   >>> for j, s in enumerate(tri.simplices):
   ...     p = points[s].mean(axis=0)
   ...     plt.text(p[0], p[1], '#%d' % j, ha='center') # label triangles
   >>> plt.xlim(-0.5, 1.5); plt.ylim(-0.5, 1.5)
   >>> plt.show()

The structure of the triangulation is encoded in the following way:
the ``simplices`` attribute contains the indices of the points in the
``points`` array that make up the triangle. For instance:

>>> i = 1
>>> tri.simplices[i,:]
array([3, 1, 0], dtype=int32)
>>> points[tri.simplices[i,:]]
array([[ 1. ,  1. ],
       [ 0. ,  1.1],
       [ 0. ,  0. ]])

Moreover, neighboring triangles can also be found out:

>>> tri.neighbors[i]
array([-1,  0, -1], dtype=int32)

What this tells us is that this triangle has triangle #0 as a neighbor,
but no other neighbors. Moreover, it tells us that neighbor 0 is
opposite the vertex 1 of the triangle:

>>> points[tri.simplices[i, 1]]
array([ 0. ,  1.1])

Indeed, from the figure we see that this is the case.

Qhull can also perform tesselations to simplices also for
higher-dimensional point sets (for instance, subdivision into
tetrahedra in 3-D).


Coplanar points
^^^^^^^^^^^^^^^

It is important to note that not *all* points necessarily appear as
vertices of the triangulation, due to numerical precision issues in
forming the triangulation.  Consider the above with a duplicated
point:

>>> points = np.array([[0, 0], [0, 1], [1, 0], [1, 1], [1, 1]])
>>> tri = Delaunay(points)
>>> np.unique(tri.simplices.ravel())
array([0, 1, 2, 3], dtype=int32)

Observe that point #4, which is a duplicate, does not occur as a
vertex of the triangulation. That this happened is recorded:

>>> tri.coplanar
array([[4, 0, 3]], dtype=int32)

This means that point 4 resides near triangle 0 and vertex 3, but is
not included in the triangulation.

Note that such degeneracies can occur not only because of duplicated
points, but also for more complicated geometrical reasons, even in
point sets that at first sight seem well-behaved.

However, Qhull has the "QJ" option, which instructs it to perturb the
input data randomly until degeneracies are resolved:

>>> tri = Delaunay(points, qhull_options="QJ Pp")
>>> points[tri.simplices]
array([[[1, 0],
        [1, 1],
        [0, 0]],
       [[1, 1],
        [1, 1],
        [1, 0]],
       [[1, 1],
        [0, 1],
        [0, 0]],
       [[0, 1],
        [1, 1],
        [1, 1]]])

Two new triangles appeared. However, we see that they are degenerate
and have zero area.


Convex hulls
------------

Convex hull is the smallest convex object containing all points in a
given point set.

These can be computed via the Qhull wrappers in `scipy.spatial` as
follows:

.. plot::

   >>> from scipy.spatial import ConvexHull
   >>> points = np.random.rand(30, 2)   # 30 random points in 2-D
   >>> hull = ConvexHull(points)
   
   The convex hull is represented as a set of N-1 dimensional simplices,
   which in 2-D means line segments. The storage scheme is exactly the
   same as for the simplices in the Delaunay triangulation discussed
   above.
   
   We can illustrate the above result:

   >>> import matplotlib.pyplot as plt
   >>> plt.plot(points[:,0], points[:,1], 'o')
   >>> for simplex in hull.simplices:
   ...     plt.plot(points[simplex,0], points[simplex,1], 'k-')
   >>> plt.show()

The same can be achieved with `scipy.spatial.convex_hull_plot_2d`.


Voronoi diagrams
----------------

A Voronoi diagram is a subdivision of the space into the nearest
neighborhoods of a given set of points.

There are two ways to approach this object using `scipy.spatial`.
First, one can use the `KDTree` to answer the question "which of the
points is closest to this one", and define the regions that way:

.. plot::

   >>> from scipy.spatial import KDTree
   >>> points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2],
   ...                    [2, 0], [2, 1], [2, 2]])
   >>> tree = KDTree(points)
   >>> tree.query([0.1, 0.1])
   (0.14142135623730953, 0)
   
   So the point ``(0.1, 0.1)`` belongs to region ``0``. In color:
   
   >>> x = np.linspace(-0.5, 2.5, 31)
   >>> y = np.linspace(-0.5, 2.5, 33)
   >>> xx, yy = np.meshgrid(x, y)
   >>> xy = np.c_[xx.ravel(), yy.ravel()]
   >>> import matplotlib.pyplot as plt
   >>> plt.pcolor(x, y, tree.query(xy)[1].reshape(33, 31))
   >>> plt.plot(points[:,0], points[:,1], 'ko')
   >>> plt.show()
   
   This does not, however, give the Voronoi diagram as a geometrical
   object.
   
   The representation in terms of lines and points can be again
   obtained via the Qhull wrappers in `scipy.spatial`:
   
   >>> from scipy.spatial import Voronoi
   >>> vor = Voronoi(points)
   >>> vor.vertices
   array([[ 0.5,  0.5],
          [ 1.5,  0.5],
          [ 0.5,  1.5],
          [ 1.5,  1.5]])
   
   The Voronoi vertices denote the set of points forming the polygonal
   edges of the Voronoi regions. In this case, there are 9 different
   regions:
   
   >>> vor.regions
   [[-1, 0], [-1, 1], [1, -1, 0], [3, -1, 2], [-1, 3], [-1, 2], [3, 1, 0, 2], [2, -1, 0], [3, -1, 1]]
   
   Negative value ``-1`` again indicates a point at infinity. Indeed,
   only one of the regions, ``[3, 1, 0, 2]``, is bounded. Note here that
   due to similar numerical precision issues as in Delaunay triangulation
   above, there may be fewer Voronoi regions than input points.
   
   The ridges (lines in 2-D) separating the regions are described as a
   similar collection of simplices as the convex hull pieces:
   
   >>> vor.ridge_vertices
   [[-1, 0], [-1, 0], [-1, 1], [-1, 1], [0, 1], [-1, 3], [-1, 2], [2, 3], [-1, 3], [-1, 2], [0, 2], [1, 3]]
   
   These numbers indicate indices of the Voronoi vertices making up the
   line segments. ``-1`` is again a point at infinity --- only four of
   the 12 lines is a bounded line segment while the others extend to
   infinity.
   
   The Voronoi ridges are perpendicular to lines drawn between the
   input points. Which two points each ridge corresponds to is also
   recorded:
   
   >>> vor.ridge_points
   array([[0, 1],
          [0, 3],
          [6, 3],
          [6, 7],
          [3, 4],
          [5, 8],
          [5, 2],
          [5, 4],
          [8, 7],
          [2, 1],
          [4, 1],
          [4, 7]], dtype=int32)
   
   This information, taken together, is enough to construct the full
   diagram.
   
   We can plot it as follows. First the points and the Voronoi vertices:
   
   >>> plt.plot(points[:, 0], points[:, 1], 'o')
   >>> plt.plot(vor.vertices[:, 0], vor.vertices[:, 1], '*')
   >>> plt.xlim(-1, 3); plt.ylim(-1, 3)
   
   Plotting the finite line segments goes as for the convex hull,
   but now we have to guard for the infinite edges:
   
   >>> for simplex in vor.ridge_vertices:
   ...     simplex = np.asarray(simplex)
   ...     if np.all(simplex >= 0):
   ...         plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-')
   
   The ridges extending to infinity require a bit more care:
   
   >>> center = points.mean(axis=0)
   >>> for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
   ...     simplex = np.asarray(simplex)
   ...     if np.any(simplex < 0):
   ...         i = simplex[simplex >= 0][0] # finite end Voronoi vertex
   ...         t = points[pointidx[1]] - points[pointidx[0]] # tangent
   ...         t /= np.linalg.norm(t)
   ...         n = np.array([-t[1], t[0]]) # normal
   ...         midpoint = points[pointidx].mean(axis=0)
   ...         far_point = vor.vertices[i] + np.sign(np.dot(midpoint - center, n)) * n * 100
   ...         plt.plot([vor.vertices[i,0], far_point[0]], 
   ...                  [vor.vertices[i,1], far_point[1]], 'k--')
   >>> plt.show()
   
This plot can also be created using `scipy.spatial.voronoi_plot_2d`.
