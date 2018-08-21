"""
Unit tests for function fitting routines.
"""

import numpy as np
from numpy.testing import (
    assert_allclose, assert_almost_equal
)
from scipy.optimize.func import (exp_fit, pow_fit)


class TestExpFit(object):
    def test_no_noise(self):
        expected, x, y = self.make_fit(331849)
        actual = exp_fit(x, y)
        assert_allclose(actual, expected, rtol=1e-3)

    def test_n_points(self):
        expected, x, y = self.make_fit(743920, n=100)
        actual_h = exp_fit(x, y)
        assert_allclose(actual_h, expected, rtol=1e-4)

        _, x, y = self.make_fit(743920, n=1000)
        actual_k = exp_fit(x, y)
        assert_allclose(actual_k, expected, rtol=1e-6)

        assert np.all(np.abs(actual_k - expected) <
                      np.abs(actual_h - expected))

    def test_noise(self):
        original, x, y = self.make_fit(789376)
        expected = (10.26939769, -0.54780143, -0.92798913)
        data = y + np.random.normal(scale=0.05 * np.ptp(y), size=x.size)
        actual = exp_fit(x, data, sorted=True)

        # Check that the fit comes out as expected
        assert_almost_equal(actual, expected)

        # Check that the RMS of the fit is better than that of the data
        fit = self.func(x, actual)
        rms_original = np.sqrt(np.square(y - data).sum())
        rms_actual = np.sqrt(np.square(fit - data).sum())
        assert rms_original > rms_actual

    def make_fit(self, seed, n=100):
        np.random.seed(seed)
        expected = (np.random.normal(scale=10),
                    np.random.normal(),
                    np.random.normal())
        x = np.linspace(-3, 7, n)
        y = self.func(x, expected)
        return expected, x, y

    @staticmethod
    def func(x, p):
        a, b, c = p
        return a + b * np.exp(c * x)


class TestPowFit(object):
    def test_no_noise(self):
        expected, x, y = self.make_fit(482293)
        actual = pow_fit(x, y, sorted=True)
        assert_allclose(actual, expected, rtol=1e-4)

    def test_n_points(self):
        expected, x, y = self.make_fit(457823, n=1000)
        actual_h = pow_fit(x, y)
        assert_allclose(actual_h, expected, rtol=1e-3)

        _, x, y = self.make_fit(457823, n=10000)
        actual_k = pow_fit(x, y)
        assert_allclose(actual_k, expected, rtol=1e-6)

        assert np.all(np.abs(actual_k - expected) <
                      np.abs(actual_h - expected))

    def test_noise(self):
        original, x, y = self.make_fit(118119)
        expected = (6.61782891, -2.53994958, -1.51472903)
        data = y + np.random.normal(scale=0.05 * np.ptp(y), size=x.size)
        actual = pow_fit(x, data, sorted=True)

        # Check that the fit comes out as expected
        assert_almost_equal(actual, expected)

        # Check that the RMS of the fit is better than that of the data
        fit = self.func(x, actual)
        rms_original = np.sqrt(np.square(y - data).sum())
        rms_actual = np.sqrt(np.square(fit - data).sum())
        assert rms_original > rms_actual

    def make_fit(self, seed, n=1000):
        np.random.seed(seed)
        expected = (np.random.normal(scale=10),
                    np.random.normal(),
                    np.random.normal())
        x = np.linspace(0.1, 10.1, n)
        y = self.func(x, expected)
        return expected, x, y

    @staticmethod
    def func(x, p):
        a, b, c = p
        return a + b * np.power(x, c)
