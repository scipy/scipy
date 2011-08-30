import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, run_module_suite

import scipy.spatial.qhull as qhull

class TestUtilities(object):
    """
    Check that utility functions work.

    """

    def test_find_simplex(self):
        # Simple check that simplex finding works
        points = np.array([(0,0), (0,1), (1,1), (1,0)], dtype=np.double)
        tri = qhull.Delaunay(points)

        # +---+
        # |\ 0|
        # | \ |
        # |1 \|
        # +---+

        assert_equal(tri.vertices, [[3, 1, 2], [3, 1, 0]])

        for p in [(0.25, 0.25, 1),
                  (0.75, 0.75, 0),
                  (0.3, 0.2, 1)]:
            i = tri.find_simplex(p[:2])
            assert_equal(i, p[2], err_msg='%r' % (p,))
            j = qhull.tsearch(tri, p[:2])
            assert_equal(i, j)

    def test_plane_distance(self):
        # Compare plane distance from hyperplane equations obtained from Qhull
        # to manually computed plane equations
        x = np.array([(0,0), (1, 1), (1, 0), (0.99189033, 0.37674127),
                      (0.99440079, 0.45182168)], dtype=np.double)
        p = np.array([0.99966555, 0.15685619], dtype=np.double)

        tri = qhull.Delaunay(x)

        z = tri.lift_points(x)
        pz = tri.lift_points(p)

        dist = tri.plane_distance(p)

        for j, v in enumerate(tri.vertices):
            x1 = z[v[0]]
            x2 = z[v[1]]
            x3 = z[v[2]]

            n = np.cross(x1 - x3, x2 - x3)
            n /= np.sqrt(np.dot(n, n))
            n *= -np.sign(n[2])

            d = np.dot(n, pz - x3)

            assert_almost_equal(dist[j], d)

    def test_convex_hull(self):
        # Simple check that the convex hull seems to works
        points = np.array([(0,0), (0,1), (1,1), (1,0)], dtype=np.double)
        tri = qhull.Delaunay(points)

        # +---+
        # |\ 0|
        # | \ |
        # |1 \|
        # +---+

        assert_equal(tri.convex_hull, [[1, 2], [3, 2], [1, 0], [3, 0]])

class TestRidgeIter2D(object):

    def _check_ridges(self, tri, vertex, expected):
        got = [(v1, v2) for v1, v2, i, t in qhull.RidgeIter2D(tri, vertex)]
        got.sort()
        expected.sort()
        assert_equal(got, expected, err_msg="%d: %r != %r" % (
            vertex, got, expected))

    def test_triangle(self):
        points = np.array([(0,0), (0,1), (1,0)], dtype=np.double)
        tri = qhull.Delaunay(points)

        # 1
        # +
        # |\
        # | \
        # |0 \
        # +---+
        # 0   2

        self._check_ridges(tri, 0, [(0, 1), (0, 2)])
        self._check_ridges(tri, 1, [(1, 0), (1, 2)])
        self._check_ridges(tri, 2, [(2, 0), (2, 1)])

    def test_rectangle(self):
        points = np.array([(0,0), (0,1), (1,1), (1,0)], dtype=np.double)
        tri = qhull.Delaunay(points)

        # 1   2
        # +---+
        # |\ 0|
        # | \ |
        # |1 \|
        # +---+
        # 0   3

        self._check_ridges(tri, 0, [(0, 1), (0, 3)])
        self._check_ridges(tri, 1, [(1, 0), (1, 3), (1, 2)])
        self._check_ridges(tri, 2, [(2, 1), (2, 3)])
        self._check_ridges(tri, 3, [(3, 0), (3, 1), (3, 2)])

    def test_complicated(self):
        points = np.array([(0,0), (0,1), (1,1), (1,0),
                           (0.5, 0.5), (0.9, 0.5)], dtype=np.double)
        tri = qhull.Delaunay(points)

        #  1                       2
        #  +-----------------------+
        #  | \-                 /-||
        #  |   \-      0      /-  /|
        #  |     \-         /-   / |
        #  |       \-     /-    |  |
        #  |         \-4/-  4  5/  |
        #  |   1       +-------+  3|
        #  |         -/  \- 5   \  |
        #  |      --/      \--   \ |
        #  |   --/     2      \- | |
        #  | -/                 \-\|
        #  +-----------------------+
        #  0                       3
        #

        self._check_ridges(tri, 0, [(0, 1), (0, 3), (0, 4)])
        self._check_ridges(tri, 1, [(1, 0), (1, 2), (1, 4)])
        self._check_ridges(tri, 2, [(2, 1), (2, 4), (2, 5), (2, 3)])
        self._check_ridges(tri, 3, [(3, 0), (3, 4), (3, 5), (3, 2)])
        self._check_ridges(tri, 4, [(4, 0), (4, 1), (4, 2), (4, 3), (4, 5)])
        self._check_ridges(tri, 5, [(5, 2), (5, 3), (5, 4)])


class TestTriangulation(object):
    """
    Check that triangulation works.

    """

    def test_nd_simplex(self):
        # simple smoke test: triangulate a n-dimensional simplex
        for nd in xrange(2, 8):
            points = np.zeros((nd+1, nd))
            for j in xrange(nd):
                points[j,j] = 1.0
            points[-1,:] = 1.0

            tri = qhull.Delaunay(points)

            tri.vertices.sort()

            assert_equal(tri.vertices, np.arange(nd+1, dtype=np.int)[None,:])
            assert_equal(tri.neighbors, -1 + np.zeros((nd+1), dtype=np.int)[None,:])

    def test_2d_square(self):
        # simple smoke test: 2d square
        points = np.array([(0,0), (0,1), (1,1), (1,0)], dtype=np.double)
        tri = qhull.Delaunay(points)

        assert_equal(tri.vertices, [[3, 1, 2], [3, 1, 0]])
        assert_equal(tri.neighbors, [[-1, -1, 1], [-1, -1, 0]])

    def test_duplicate_points(self):
        x  = np.array([0, 1, 0, 1], dtype=np.float64)
        y  = np.array([0, 0, 1, 1], dtype=np.float64)

        xp = np.r_[x, x]
        yp = np.r_[y, y]

        # shouldn't fail on duplicate points
        tri = qhull.Delaunay(np.c_[x, y])
        tri2 = qhull.Delaunay(np.c_[xp, yp])

    pathological_data_1 = np.array([
        [-3.14,-3.14], [-3.14,-2.36], [-3.14,-1.57], [-3.14,-0.79],
        [-3.14,0.0], [-3.14,0.79], [-3.14,1.57], [-3.14,2.36],
        [-3.14,3.14], [-2.36,-3.14], [-2.36,-2.36], [-2.36,-1.57],
        [-2.36,-0.79], [-2.36,0.0], [-2.36,0.79], [-2.36,1.57],
        [-2.36,2.36], [-2.36,3.14], [-1.57,-0.79], [-1.57,0.79],
        [-1.57,-1.57], [-1.57,0.0], [-1.57,1.57], [-1.57,-3.14],
        [-1.57,-2.36], [-1.57,2.36], [-1.57,3.14], [-0.79,-1.57],
        [-0.79,1.57], [-0.79,-3.14], [-0.79,-2.36], [-0.79,-0.79],
        [-0.79,0.0], [-0.79,0.79], [-0.79,2.36], [-0.79,3.14],
        [0.0,-3.14], [0.0,-2.36], [0.0,-1.57], [0.0,-0.79], [0.0,0.0],
        [0.0,0.79], [0.0,1.57], [0.0,2.36], [0.0,3.14], [0.79,-3.14],
        [0.79,-2.36], [0.79,-0.79], [0.79,0.0], [0.79,0.79],
        [0.79,2.36], [0.79,3.14], [0.79,-1.57], [0.79,1.57],
        [1.57,-3.14], [1.57,-2.36], [1.57,2.36], [1.57,3.14],
        [1.57,-1.57], [1.57,0.0], [1.57,1.57], [1.57,-0.79],
        [1.57,0.79], [2.36,-3.14], [2.36,-2.36], [2.36,-1.57],
        [2.36,-0.79], [2.36,0.0], [2.36,0.79], [2.36,1.57],
        [2.36,2.36], [2.36,3.14], [3.14,-3.14], [3.14,-2.36],
        [3.14,-1.57], [3.14,-0.79], [3.14,0.0], [3.14,0.79],
        [3.14,1.57], [3.14,2.36], [3.14,3.14],
    ])

    pathological_data_2 = np.array([
        [-1, -1                          ], [-1, 0], [-1, 1],
        [ 0, -1                          ], [ 0, 0], [ 0, 1],
        [ 1, -1 - np.finfo(np.float_).eps], [ 1, 0], [ 1, 1],
    ])

    def test_pathological(self):
        # both should succeed
        tri = qhull.Delaunay(self.pathological_data_1)
        assert_equal(tri.points[tri.vertices].max(),
                     self.pathological_data_1.max())
        assert_equal(tri.points[tri.vertices].min(),
                     self.pathological_data_1.min())

        tri = qhull.Delaunay(self.pathological_data_2)
        assert_equal(tri.points[tri.vertices].max(),
                     self.pathological_data_2.max())
        assert_equal(tri.points[tri.vertices].min(),
                     self.pathological_data_2.min())

if __name__ == "__main__":
    run_module_suite()
