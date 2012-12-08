import os
import sys

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, run_module_suite,\
     assert_, dec, assert_allclose, assert_array_equal

import scipy.spatial.qhull as qhull
from scipy.spatial import cKDTree as KDTree

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

    def _check_barycentric_transforms(self, tri, err_msg="",
                                      unit_cube=False,
                                      unit_cube_tol=0):
        """Check that a triangulation has reasonable barycentric transforms"""
        vertices = tri.points[tri.vertices]
        sc = 1/(tri.ndim + 1.0)
        centroids = vertices.sum(axis=1) * sc

        # Either: (i) the simplex has a `nan` barycentric transform,
        # or, (ii) the centroid is in the simplex

        def barycentric_transform(tr, x):
            ndim = tr.shape[1]
            r = tr[:,-1,:]
            Tinv = tr[:,:-1,:]
            return np.einsum('ijk,ik->ij', Tinv, x - r)

        eps = np.finfo(float).eps

        c = barycentric_transform(tri.transform, centroids)
        olderr = np.seterr(invalid="ignore")
        try:
            ok = np.isnan(c).all(axis=1) | (abs(c - sc)/sc < 0.1).all(axis=1)
        finally:
            np.seterr(**olderr)

        assert_(ok.all(), "%s %s" % (err_msg, np.where(~ok)))

        # Invalid simplices must be (nearly) zero volume
        q = vertices[:,:-1,:] - vertices[:,-1,None,:]
        volume = np.array([np.linalg.det(q[k,:,:])
                           for k in range(tri.nsimplex)])
        ok = np.isfinite(tri.transform[:,0,0]) | (volume < np.sqrt(eps))
        assert_(ok.all(), "%s %s" % (err_msg, np.where(~ok)))

        # Also, find_simplex for the centroid should end up in some
        # simplex for the non-degenerate cases
        j = tri.find_simplex(centroids)
        ok = (j != -1) | np.isnan(tri.transform[:,0,0])
        assert_(ok.all(), "%s %s" % (err_msg, np.where(~ok)))

        if unit_cube:
            # If in unit cube, no interior point should be marked out of hull
            at_boundary =  (centroids <= unit_cube_tol).any(axis=1)
            at_boundary |= (centroids >= 1 - unit_cube_tol).any(axis=1)

            ok = (j != -1) | at_boundary
            assert_(ok.all(), "%s %s" % (err_msg, np.where(~ok)))

    @dec.skipif(np.version.short_version < '1.6', "No einsum in numpy 1.5.x")
    def test_degenerate_barycentric_transforms(self):
        # The triangulation should not produce invalid barycentric
        # transforms that stump the simplex finding
        data = np.load(os.path.join(os.path.dirname(__file__), 'data',
                                    'degenerate_pointset.npz'))
        points = data['c']
        data.close()

        tri = qhull.Delaunay(points)

        # Check that there are not too many invalid simplices
        bad_count = np.isnan(tri.transform[:,0,0]).sum()
        assert_(bad_count < 20, bad_count)

        # Check the transforms
        self._check_barycentric_transforms(tri)

    @dec.slow
    @dec.skipif(np.version.short_version < '1.6', "No einsum in numpy 1.5.x")
    def test_more_barycentric_transforms(self):
        # Triangulate some "nasty" grids

        eps = np.finfo(float).eps

        npoints = {2: 70, 3: 11, 4: 5, 5: 3}

        for ndim in xrange(2, 6):
            # Generate an uniform grid in n-d unit cube
            x = np.linspace(0, 1, npoints[ndim])
            grid = np.c_[map(np.ravel, np.broadcast_arrays(*np.ix_(*([x]*ndim))))].T

            err_msg = "ndim=%d" % ndim

            # Check using regular grid
            tri = qhull.Delaunay(grid)
            self._check_barycentric_transforms(tri, err_msg=err_msg,
                                               unit_cube=True)

            # Check with eps-perturbations
            np.random.seed(1234)
            m = (np.random.rand(grid.shape[0]) < 0.2)
            grid[m,:] += 2*eps*(np.random.rand(*grid[m,:].shape) - 0.5)

            tri = qhull.Delaunay(grid)
            self._check_barycentric_transforms(tri, err_msg=err_msg,
                                               unit_cube=True,
                                               unit_cube_tol=2*eps)

            # Check with duplicated data
            tri = qhull.Delaunay(np.r_[grid, grid])
            self._check_barycentric_transforms(tri, err_msg=err_msg,
                                               unit_cube=True,
                                               unit_cube_tol=2*eps)

            # Check with larger perturbations
            np.random.seed(4321)
            m = (np.random.rand(grid.shape[0]) < 0.2)
            grid[m,:] += 1000*eps*(np.random.rand(*grid[m,:].shape) - 0.5)

            tri = qhull.Delaunay(grid)
            self._check_barycentric_transforms(tri, err_msg=err_msg,
                                               unit_cube=True,
                                               unit_cube_tol=1500*eps)

            # Check with yet larger perturbations
            np.random.seed(4321)
            m = (np.random.rand(grid.shape[0]) < 0.2)
            grid[m,:] += 1e6*eps*(np.random.rand(*grid[m,:].shape) - 0.5)

            tri = qhull.Delaunay(grid)
            self._check_barycentric_transforms(tri, err_msg=err_msg,
                                               unit_cube=True,
                                               unit_cube_tol=1.5e6*eps)


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


class TestDelaunay(object):
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

    def test_joggle(self):
        # Check that the option QJ indeed guarantees that all input points
        # occur as vertices of the triangulation

        points = np.random.rand(10, 2)
        points = np.r_[points, points] # duplicate input data

        tri = qhull.Delaunay(points, qhull_options="QJ Qbb Pp")
        assert_array_equal(np.unique(tri.simplices.ravel()),
                           np.arange(len(points)))

    def test_coplanar(self):
        # Check that the coplanar point output option indeed works
        points = np.random.rand(10, 2)
        points = np.r_[points, points] # duplicate input data

        tri = qhull.Delaunay(points)

        assert_(len(np.unique(tri.simplices.ravel())) == len(points)//2)
        assert_(len(tri.coplanar) == len(points)//2)

        assert_(len(np.unique(tri.coplanar[:,2])) == len(points)//2)

        assert_(np.all(tri.vertex_to_simplex >= 0))


class TestConvexHull:
    def test_hull_consistency_tri(self):
        # Check that a convex hull returned by qhull in ndim
        # and the hull constructed from ndim delaunay agree
        np.random.seed(1234)

        datasets = {'pathological_1': TestDelaunay.pathological_data_1,
                    'pathological_2': TestDelaunay.pathological_data_2}
        for nd in range(2, 8):
            points = np.random.rand(30, nd)
            datasets['random-%dd' % nd] = points

        def check(name):
            points = datasets[name]

            tri = qhull.Delaunay(points)
            hull = qhull.ConvexHull(points)

            def sorted_tuple(x):
                return tuple(sorted(x))

            facets_1 = set(map(sorted_tuple, tri.convex_hull.tolist()))
            facets_2 = set(map(sorted_tuple, hull.simplices.tolist()))

            if facets_1 == facets_2:
                return # OK

            if points.shape[1] == 2:
                # The direct check fails for the pathological cases
                # --- then the convex hull from Delaunay differs (due
                # to rounding error etc.) from the hull computed
                # otherwise, by the question whether (tricoplanar)
                # points that lie almost exactly on the hull are
                # included as vertices of the hull or not.
                #
                # So we check the result, and accept it if the Delaunay
                # hull line segments are a subset of the usual hull.

                eps = 1000 * np.finfo(float).eps

                for a, b in facets_1:
                    for ap, bp in facets_2:
                        t = points[bp] - points[ap]
                        t /= np.linalg.norm(t)       # tangent
                        n = np.array([-t[1], t[0]])  # normal

                        # check that the two line segments are parallel
                        # to the same line
                        c1 = np.dot(n, points[b] - points[ap])
                        c2 = np.dot(n, points[a] - points[ap])
                        if not np.allclose(np.dot(c1, n), 0):
                            continue
                        if not np.allclose(np.dot(c2, n), 0):
                            continue

                        # Check that the segment (a, b) is contained in (ap, bp)
                        c1 = np.dot(t, points[a] - points[ap])
                        c2 = np.dot(t, points[b] - points[ap])
                        c3 = np.dot(t, points[bp] - points[ap])
                        if c1 < -eps or c1 > c3 + eps:
                            continue
                        if c2 < -eps or c2 > c3 + eps:
                            continue

                        # OK:
                        break
                    else:
                        raise AssertionError("comparison fails")
                return
            raise AssertionError("comparison fails...")

        for name in datasets.keys():
            yield check, name

class TestVoronoi:
    def test_simple(self):
        # Simple case with known Voronoi diagram
        points = [(0, 0), (0, 1), (0, 2),
                  (1, 0), (1, 1), (1, 2),
                  (2, 0), (2, 1), (2, 2)]

        # qhull v o Fv Qbb Qc Qz < dat
        output = """
        2
        5 10 1
        -10.101 -10.101 
           0.5    0.5 
           1.5    0.5 
           0.5    1.5 
           1.5    1.5 
        2 0 1
        3 3 0 1
        2 0 3
        3 2 0 1
        4 4 3 1 2
        3 4 0 3
        2 0 2
        3 4 0 2
        2 0 4
        0
        12
        4 0 3 0 1
        4 0 1 0 1
        4 1 4 1 3
        4 1 2 0 3
        4 2 5 0 3
        4 3 4 1 2
        4 3 6 0 2
        4 4 5 3 4
        4 4 7 2 4
        4 5 8 0 4
        4 6 7 0 2
        4 7 8 0 4
        """
        self._compare_qvoronoi(points, output)

    def _compare_qvoronoi(self, points, output):
        """Compare to output from 'qvoronoi o Fv < data' to Voronoi()"""

        # Parse output
        output = map(lambda x: map(float, x.split()),
                     output.strip().splitlines())
        nvertex = int(output[1][0])
        vertices = map(tuple, output[3:2+nvertex]) # exclude inf
        nregion = int(output[1][1])
        regions = [[int(y)-1 for y in x[1:]]
                   for x in output[2+nvertex:2+nvertex+nregion]]
        nridge = int(output[2+nvertex+nregion][0])
        ridge_points = [[int(y) for y in x[1:3]]
                        for x in output[3+nvertex+nregion:]]
        ridge_vertices = [[int(y)-1 for y in x[3:]]
                          for x in output[3+nvertex+nregion:]]

        # Compare results
        vor = qhull.Voronoi(points)

        def sorttuple(x):
            return tuple(sorted(x))

        assert_allclose(vor.vertices, vertices)
        assert_equal(set(map(tuple, vor.regions)),
                     set(map(tuple, regions)))

        p1 = zip(map(sorttuple, ridge_points), map(sorttuple, ridge_vertices))
        p2 = zip(map(sorttuple, vor.ridge_points.tolist()),
                 map(sorttuple, vor.ridge_vertices))
        p1.sort()
        p2.sort()

        print "\n".join(map(repr, p1))
        print "--"
        print "\n".join(map(repr, p2))

        assert_equal(p1, p2)
    
    def test_ridges(self):
        # Check that the ridges computed by Voronoi indeed separate
        # the regions of nearest neighborhood, by comparing the result
        # to KDTree.

        np.random.seed(1234)

        datasets = {'pathological_1': TestDelaunay.pathological_data_1,
                    'pathological_2': TestDelaunay.pathological_data_2}
        for nd in range(2, 8):
            points = np.random.rand(30, nd)
            datasets['random-%dd' % nd] = points

        def check(name):
            points = datasets[name]

            tree = KDTree(points)
            vor = qhull.Voronoi(points)

            for p, v in vor.ridge_dict.items():
                # consider only finite ridges
                if not np.all(np.asarray(v) >= 0):
                    continue

                ridge_midpoint = vor.vertices[v].mean(axis=0)
                d = 1e-6 * (points[p[0]] - ridge_midpoint)

                dist, k = tree.query(ridge_midpoint + d, k=1)
                assert_equal(k, p[0])

                dist, k = tree.query(ridge_midpoint - d, k=1)
                assert_equal(k, p[1])

        for name in datasets.keys():
            yield check, name

if __name__ == "__main__":
    run_module_suite()
