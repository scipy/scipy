import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial import AlphaShapes, ConvexHull


TOL = 1E-10


class TestAlphaShapes(object):

    def setup_method(self):
        seed = 0
        size = 100
        dim = 3
        rng = np.random.RandomState(seed=seed)
        self.points = rng.normal(size=(size, dim))
        self.alpha = AlphaShapes(self.points)

    def test_circumcenters(self):
        # test vertex distances from circumcenters are equal
        alpha = self.alpha
        tetrahedra = alpha.points[alpha.simplices]
        for t, c in zip(tetrahedra, alpha.circumcenters):
            delta = t - c
            sizes = np.linalg.norm(delta, axis=1)
            assert_allclose(sizes, sizes[0], atol=TOL)

    def test_radii(self):
        alpha = self.alpha
        deltas = alpha.circumcenters - alpha.points[alpha.simplices[:, 0]]
        radii = np.linalg.norm(deltas, axis=1)
        assert_allclose(radii, alpha.radii, atol=TOL)

    def test_boundary_facet_intervals(self):
        alpha = self.alpha
        assert (alpha._start <= alpha._end).all()

    def test_convex_hull(self):
        alpha = self.alpha
        simplices = np.sort(ConvexHull(self.points).simplices)
        facets = np.sort(alpha.get_boundary_facets(float("inf")))
        assert simplices.shape == facets.shape

        findices = np.lexsort(facets.T)
        sindices = np.lexsort(simplices.T)
        assert (facets[findices] == simplices[sindices]).all()

    def test_zero_boundary_radius(self):
        alpha = self.alpha
        assert alpha.get_boundary_facets(0).shape == (0,)

    def test_connected_components(self):
        alpha = self.alpha
        threshold = alpha.thresholds[1]
        indices = np.where(alpha.radii <= threshold)[0]
        hits = np.unique(alpha.simplices[indices])
        assert (hits == np.arange(len(self.points))).all()
