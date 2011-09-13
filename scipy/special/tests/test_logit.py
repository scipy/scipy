import numpy as np
from numpy.testing import TestCase, assert_equal, assert_almost_equal
from scipy.special import logit, expit


class TestLogit(TestCase):
    def check_logit_out(self, dtype, expected):
        a = np.linspace(0,1,10)
        a = np.array(a, dtype=dtype)
        olderr = np.seterr(divide='ignore')
        try:
            actual = logit(a)
        finally:
            np.seterr(**olderr)

        if np.__version__ >= '1.6':
            assert_almost_equal(actual, expected)
        else:
            assert_almost_equal(actual[1:-1], expected[1:-1])

        assert_equal(actual.dtype, np.dtype(dtype))

    def test_float32(self):
        expected = np.array([-np.inf, -2.07944155,
                            -1.25276291, -0.69314718,
                            -0.22314353, 0.22314365,
                            0.6931473 ,  1.25276303,
                            2.07944155, np.inf], dtype=np.float32)
        self.check_logit_out('f4', expected)

    def test_float64(self):
        expected = np.array([-np.inf, -2.07944154,
                            -1.25276297, -0.69314718,
                            -0.22314355, 0.22314355,
                            0.69314718,  1.25276297,
                            2.07944154,         np.inf])
        self.check_logit_out('f8', expected)

    def test_nan(self):
        expected = np.array([np.nan]*4)
        olderr = np.seterr(invalid='ignore')
        try:
            actual = logit(np.array([-3., -2., 2., 3.]))
        finally:
            np.seterr(**olderr)

        assert_equal(expected, actual)


class TestExpit(TestCase):
    def check_expit_out(self, dtype, expected):
        a = np.linspace(-4,4,10)
        a = np.array(a, dtype=dtype)
        actual = expit(a)
        assert_almost_equal(actual, expected)
        assert_equal(actual.dtype, np.dtype(dtype))

    def test_float32(self):
        expected = np.array([ 0.01798621,  0.04265125,
                            0.09777259,  0.20860852,
                            0.39068246, 0.60931754,
                            0.79139149,  0.9022274 ,
                            0.95734876,  0.98201376], dtype=np.float32)
        self.check_expit_out('f4',expected)

    def test_float64(self):
        expected = np.array([ 0.01798621,  0.04265125,
                            0.0977726 ,  0.20860853,
                            0.39068246, 0.60931754,
                            0.79139147,  0.9022274 ,
                            0.95734875,  0.98201379])
        self.check_expit_out('f8', expected)


