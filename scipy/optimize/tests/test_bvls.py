from scipy.optimize import bounded_lstsq
from numpy import arange, dot, float32
from numpy.linalg import norm
from numpy.testing import (TestCase, assert_array_almost_equal, assert_equal,
                           assert_)

class TestBVLS(TestCase):

    def test_bvls(self):
        X = arange(25.0, dtype=float32).reshape(-1, 5)
        b = arange(5.0, dtype=float32)
        y = dot(X, b)
        bb, res, nsetp = bounded_lstsq(X, y, bounds=[(0, None), (1, None),
                                                     (2, None), (3, None),
                                                     (4, None)])
        assert_equal(b, bb)
        assert_(res < 1e-7)
        assert_(norm(dot(X,bb)-y) < 1e-7)

        bb, res, nsetp = bounded_lstsq(X, y, bounds=[(None, None),
                                                     (None, None),
                                                     (None, None),
                                                     (1, 2),
                                                     (2, 3)])
        assert_(bb[3] == 1)
        assert_(bb[4] == 3)

        # smoke test
        bb, res, nsetp = bounded_lstsq(X, y, bounds=None)
