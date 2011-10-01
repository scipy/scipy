import unittest
from numpy import array, arange, exp, linspace, empty, sqrt
from numpy.random import standard_normal
from numpy.testing import assert_equal, assert_almost_equal
from scipy.optimize import Fit, InvalidParameter

class TestGaussian(unittest.TestCase, Fit):
    def setUp(self):
        self.complain = False

    def func(self, p):
        if p[2] < 0 and self.complain:
            raise InvalidParameter(2)
        return p[0] + p[1] * exp(-((self.x - p[3]) / p[2]) ** 2 / 2) - self.y

    def test_gaussian(self):
        self.y = array([0.1, 1.5, 2.7, 3.2, 4.8, 5.3, 4.4, 3.6, 2.9, 1.0])
        self.x = arange(10.)
        r = self.fit([1., 2, 3, 4])
        assert_almost_equal(r, [-1.91, 6.87, 3.22, 4.99], decimal=2)

    def test_random(self):
        self.x = linspace(-20, 20, 100)
        self.y = standard_normal(len(self.x)) / 10
        p = array([3., 5., 3., 4.])
        self.y = self.func(p)
        pp = self.fit(p)
        assert_almost_equal(p, pp, decimal=1)

    def test_error(self):
        self.x = linspace(-20, 20, 100)
        p = array([3., 5., 3., 4.])
        err = 0.03
        M = 100
        pp = empty((M, len(p)))
        for i in range(M):
            self.y = standard_normal(len(self.x)) * err
            self.y = self.func(p)
            pp[i, :] = self.fit(p)
        assert_almost_equal(pp.std(axis=0), self.error() * err, decimal=1)
        assert_almost_equal(sqrt((self.func(p) ** 2).mean()), err, decimal=1)

    def test_negative(self):
        self.x = linspace(-20, 20, 100)
        self.y = standard_normal(len(self.x)) / 100
        p = array([1., 2., 3., 4.])
        self.y = self.func(p)
        p = array([1., 1., 20., 1.])
        self.complain = False
        p1 = self.fit(p)
        self.complain = True
        p2 = self.fit(p)
        self.complain = False
        p2[2] *= -1
        assert_almost_equal(p1, p2, decimal=5)

    def test_dependent(self):
        self.x = linspace(-20, 20, 100)
        self.y = standard_normal(len(self.x)) / 100
        p = array([1., 2., 3., 4.])
        self.y = self.func(p)
        p = array([1., 1., 20., 1., 0])
        self.complain = True
        pp = self.fit(p)
        self.complain = False
        assert_almost_equal(pp, [1, 2, 3, 4, 0], decimal=2)
        assert_almost_equal(self.error(), [0.10, 0.33, 0.66, 0.59, 0],
            decimal=1)
        

if __name__ == '__main__':
    unittest.main()
