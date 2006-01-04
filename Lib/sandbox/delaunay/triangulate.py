# 2.3 compatibility
try:
    set
except NameError:
    import sets
    set = sets.Set

import numpy as sp

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

      hull -- list of point_id's giving the nodes which form the convex hull
        of the point set. This list is sorted in counter-clockwise order.
    """
    def __init__(self, x, y):
        self.x = sp.asarray(x, dtype=sp.float64)
        self.y = sp.asarray(y, dtype=sp.float64)

        if self.x.shape != self.y.shape or len(self.x.shape) != 1:
            raise ValueError("x,y must be equal-length 1-D arrays")

        self.circumcenters, self.edge_db, self.triangle_nodes, \
            self.triangle_neighbors = delaunay(self.x, self.y)

        self.hull = self._compute_convex_hull()

    def _compute_convex_hull(self):
        """Extract the convex hull from the triangulation information.

        The output will be a list of point_id's in counter-clockwise order
        forming the convex hull of the data set.
        """
        border = (self.triangle_neighbors == -1)

        edges = {}
        edges.update(zip(self.triangle_nodes[border[:,0]][:,1],
                         self.triangle_nodes[border[:,0]][:,2]))
        edges.update(zip(self.triangle_nodes[border[:,1]][:,2],
                         self.triangle_nodes[border[:,1]][:,0]))
        edges.update(zip(self.triangle_nodes[border[:,2]][:,0],
                         self.triangle_nodes[border[:,2]][:,1]))
        
        # Take an arbitrary starting point and its subsequent node
        hull = list(edges.popitem())
        while edges:
            hull.append(edges.pop(hull[-1]))

        # hull[-1] == hull[0], so remove hull[-1]
        hull.pop()

        return hull

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

    def prep_extrapolator(self, z, bbox=None):
        if bbox is None:
            bbox = (self.x[0], self.x[0], self.y[0], self.y[0])
        minx, maxx, miny, maxy = sp.asarray(bbox, sp.float64)
        minx = min(minx, sp.minimum.reduce(self.x))
        miny = min(miny, sp.minimum.reduce(self.y))
        maxx = max(maxx, sp.maximum.reduce(self.x))
        maxy = max(maxy, sp.maximum.reduce(self.y))
        M = max((maxx-minx)/2, (maxy-miny)/2)
        midx = (minx + maxx)/2.0
        midy = (miny + maxy)/2.0

        xp, yp= sp.array([[midx+3*M, midx, midx-3*M],
                          [midy, midy+3*M, midy-3*M]])
        x1 = sp.hstack((self.x, xp))
        y1 = sp.hstack((self.y, yp))
        newtri = self.__class__(x1, y1)

        # do a least-squares fit to a plane to make pseudo-data
        xy1 = sp.ones((len(self.x), 3), sp.float64)
        xy1[:,0] = self.x
        xy1[:,1] = self.y
        from scipy import linalg
        c, res, rank, s = linalg.lstsq(xy1, z)
        zp = sp.hstack((z, xp*c[0] + yp*c[1] + c[2]))

        return newtri, zp

    def nn_extrapolator(self, z, bbox=None, default_value=sp.nan):
        newtri, zp = self.prep_extrapolator(z, bbox)
        return newtri.nn_interpolator(zp, default_value)

    def linear_extrapolator(self, z, bbox=None, default_value=sp.nan):
        newtri, zp = self.prep_extrapolator(z, bbox)
        return newtri.linear_interpolator(zp, default_value)

    def node_graph(self):
        """Return a graph of node_id's pointing to node_id's.

        The arcs of the graph correspond to the edges in the triangulation.

        {node_id: set([node_id, ...]), ...}
        """
        g = {}
        for i, j in self.edge_db:
            s = g.setdefault(i, set())
            s.add(j)
            s = g.setdefault(j, set())
            s.add(i)
        return g


