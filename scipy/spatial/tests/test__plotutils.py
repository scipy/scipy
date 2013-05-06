from __future__ import division, print_function, absolute_import

from numpy.testing import dec, assert_, assert_array_equal

try:
    import matplotlib
    matplotlib.rcParams['backend'] = 'Agg'
    import matplotlib.pyplot as plt
    has_matplotlib = True
except:
    has_matplotlib = False

from scipy.spatial import \
     delaunay_plot_2d, voronoi_plot_2d, convex_hull_plot_2d, \
     Delaunay, Voronoi, ConvexHull


class TestPlotting:
    points = [(0,0), (0,1), (1,0), (1,1)]

    @dec.skipif(not has_matplotlib, "Matplotlib not available")
    def test_delaunay(self):
        # Smoke test
        fig = plt.figure()
        obj = Delaunay(self.points)
        s_before = obj.simplices.copy()
        r = delaunay_plot_2d(obj, ax=fig.gca())
        assert_array_equal(obj.simplices, s_before)  # shouldn't modify
        assert_(r is fig)
        delaunay_plot_2d(obj, ax=fig.gca())

    @dec.skipif(not has_matplotlib, "Matplotlib not available")
    def test_voronoi(self):
        # Smoke test
        fig = plt.figure()
        obj = Voronoi(self.points)
        r = voronoi_plot_2d(obj, ax=fig.gca())
        assert_(r is fig)
        voronoi_plot_2d(obj)

    @dec.skipif(not has_matplotlib, "Matplotlib not available")
    def test_convex_hull(self):
        # Smoke test
        fig = plt.figure()
        tri = ConvexHull(self.points)
        r = convex_hull_plot_2d(tri, ax=fig.gca())
        assert_(r is fig)
        convex_hull_plot_2d(tri)
