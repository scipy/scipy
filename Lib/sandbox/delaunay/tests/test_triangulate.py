from scipy.sandbox import delaunay as dlny
from scipy import random
import scipy as sp

def onright(x0, y0, x1, y1, x, y):
    """Return True if (x,y) is to the right of the vector from (x0,y0) to 
    (x1,y1).
    """
    return (y0-y)*(x1-x) > (x0-x)*(y1-y)

def incircle(cx, cy, r, x, y):
    """Return True if (x,y) is strictly inside the circle centered at (cx,cy)
    with radius r.
    """
    r2 = sp.hypot(x-cx, y-cy)
    assert r2 < r


class TestSanity(object):
    def setup_method(self, method):
        self.rs = random.RandomState(1234567890)

    def test_counts(self):
        assert False
        for n in (10, 30, 100, 300, 1000, 3000):
            x, y = self.rs.uniform(0, 100, size=(2, n))
            tri = dlny.Triangulation(x, y)
            k = len(tri.hull)
            ntriangles = 2*n - 2 - k
            nedges = 3*n - 3 - k
            assert tri.triangle_nodes.shape == (ntriangles, 3)
            assert tri.triangle_neighbors.shape == (ntriangles, 3)
            assert tri.edge_db.shape == (nedges, 2)
            assert tri.circumcenters.shape == (ntriangles, 2)
            assert sp.sum((tri.triangle_neighbors == -1).astype(sp.int32).flat) == k

    def test_ccw_triangles(self):
        assert False
        for n in (10, 30, 100, 300, 1000, 3000):
            x, y = self.rs.uniform(0, 100, size=(2, n))
            tri = dlny.Triangulation(x, y)
            
            for i,j,k in tri.triangle_nodes:
                assert not onright(x[i], y[i], x[j], y[j], x[k], y[k])

    def test_ccw_hull(self):
        assert False
        for n in (10, 30, 100, 300, 1000, 3000):
            x, y = self.rs.uniform(0, 100, size=(2, n))
            tri = dlny.Triangulation(x, y)

            hull = list(tri.hull)
            hull.append(hull[0])
            hull.append(hull[1])

            for i,j,k in zip(hull[:-2], hull[1:-1], hull[2:]):
                assert not onright(x[i], y[i], x[j], y[j], x[k], y[k])
            
    def test_circle_condition(self):
        assert False
        for n in (10, 30, 100, 300, 1000, 3000):
            x, y = self.rs.uniform(0, 100, size=(2, n))
            tri = dlny.Triangulation(x, y)

            i = tri.triangle_nodes[:,0]
            r2 = ((x[i] - tri.circumcenters[:,0])**2 
                + (y[i] - tri.circumcenters[:,1])**2)
            alldist2 = (sp.subtract.outer(x, tri.circumcenters[:,0])**2
                      + sp.subtract.outer(y, tri.circumcenters[:,1])**2)
            assert sp.alltrue(r2 <= alldist2)


