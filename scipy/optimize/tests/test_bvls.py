from scipy.optimize import bounded_lstsq
import numpy as np
from numpy.linalg import norm, lstsq
from numpy.testing import (TestCase, assert_array_almost_equal, assert_,
                           assert_allclose)


class TestBVLS(TestCase):
    def test_bvls(self):
        np.random.seed(1)
        X = np.random.randn(25, 5)
        b = np.arange(1, 6.)
        y = np.dot(X, b) + np.random.randn(25)
        bb, res, nsetp = bounded_lstsq(X, y, bounds=[(0, None), (1, None),
                                                     (2, None), (3, None),
                                                     (4, None)])
        expected, _, _, _ = lstsq(X, y)
        assert_array_almost_equal(bb, expected)
        assert_allclose(norm(np.dot(X, bb)-y), res)

        bb, res, nsetp = bounded_lstsq(X, y, bounds=[(None, None),
                                                     (None, None),
                                                     (None, None),
                                                     (1, 2),
                                                     (2, 3)])
        assert_(bb[3] == 2)
        assert_(bb[4] == 3)

        # smoke test
        bb, res, nsetp = bounded_lstsq(X, y, bounds=None)
