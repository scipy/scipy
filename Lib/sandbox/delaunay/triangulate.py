import scipy.base as sp

from _delaunay import delaunay
from interpolate import LinearInterpolator, NNInterpolator

__all__ = ['Triangulation']

class Triangulation(object):
    """A Delaunay triangulation of points in a plane.

    Triangulation(x, y)
    x, y -- the coordinates of the points as 1-D arrays of floats

    Let us make the following definitions:
        npoints = number of points input
        nedges = number of edges in the triangulation
        ntriangles = number of triangles in the triangulation

        point_id = an integer identifying a particular point (specifically, an
            index into x and y), range(0, npoints)
        edge_id = an integer identifying a particular edge, range(0, nedges)
        triangle_id = an integer identifying a particular triangle 
            range(0, ntriangles)

    Attributes: (all should be treated as read-only to maintain consistency)
      x, y -- the coordinates of the points as 1-D arrays of floats.

      circumcenters -- (ntriangles, 2) array of floats giving the (x,y) 
        coordinates of the circumcenters of each triangle (indexed by a 
        triangle_id).

      edge_db -- (nedges, 2) array of point_id's giving the points forming 
        each edge in no particular order; indexed by an edge_id.

      triangle_nodes -- (ntriangles, 3) array of point_id's giving the points 
        forming each triangle in counter-clockwise order; indexed by a 
        triangle_id.

      triangle_neighbors -- (ntriangles, 3) array of triangle_id's giving the 
        neighboring triangle; indexed by a triangle_id. 
        
        The value can also be -1 meaning that that edge is on the convex hull of
        the points and there is no neighbor on that edge. The values are ordered
        such that triangle_neighbors[tri, i] corresponds with the edge
        *opposite* triangle_nodes[tri, i]. As such, these neighbors are also in
        counter-clockwise order.

    """
    def __init__(self, x, y):
        self.x = sp.asarray(x, dtype=sp.float64)
        self.y = sp.asarray(y, dtype=sp.float64)

        if self.x.shape != self.y.shape or len(self.x.shape) != 1:
            raise ValueError("x,y must be equal-length 1-D arrays")

        self.circumcenters, self.edge_db, self.triangle_nodes, \
            self.triangle_neighbors = delaunay(self.x, self.y)

    def linear_interpolator(self, z, default_value=sp.nan):
        """Get an object which can interpolate within the convex hull by 
        assigning a plane to each triangle.

        z -- an array of floats giving the known function values at each point 
          in the triangulation.
        """
        z = sp.asarray(z, dtype=sp.float64)
        if z.shape != self.x.shape:
            raise ValueError("z must be the same shape as x and y")

        return LinearInterpolator(self, z, default_value)

    def nn_interpolator(self, z, default_value=sp.nan):
        """Get an object which can interpolate within the convex hull by 
        the natural neighbors method.

        z -- an array of floats giving the known function values at each point 
          in the triangulation.
        """
        z = sp.asarray(z, dtype=sp.float64)
        if z.shape != self.x.shape:
            raise ValueError("z must be the same shape as x and y")

        return NNInterpolator(self, z, default_value)

