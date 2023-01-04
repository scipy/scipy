import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
import pytest

from scipy.interpolate._pava_pybind import pava
from scipy.interpolate import isotonic_regression

class TestIsotonicRegression:
    def test_simple_pava(self):
        # Test case of Busing 2020
        # https://doi.org/10.18637/jss.v102.c01
        y = np.array([8, 4, 8, 2, 2, 0, 8], dtype=np.float64)
        w = np.ones_like(y)
        r = np.full(shape=y.shape[0] + 1, fill_value=-1, dtype=np.intp)
        pava(y, w, r)
        assert_almost_equal(y, [4, 4, 4, 4, 4, 4, 8])
        # Only first 2 elements of w are changed.
        assert_almost_equal(w, [6, 1, 1, 1, 1, 1, 1])
        # Only first 3 elements of r are changed.
        assert_almost_equal(r, [0, 6, 7, -1, -1, -1, -1, -1])

    @pytest.mark.parametrize("w", [None, np.ones(7)])
    def test_simple_isotonic_regression(self, w):
        # Test case of Busing 2020
        # https://doi.org/10.18637/jss.v102.c01
        y = np.array([8, 4, 8, 2, 2, 0, 8], dtype=np.float64)
        x, w, r = isotonic_regression(y, w)
        assert_almost_equal(x, [4, 4, 4, 4, 4, 4, 8])
        assert_almost_equal(w, [6, 1])
        assert_almost_equal(r, [0, 6, 7])
        # Assert that y was not overwritten
        assert_equal(y, np.array([8, 4, 8, 2, 2, 0, 8], dtype=np.float64))
    
    @pytest.mark.parametrize("increasing", [True, False])
    def test_linspace(self, increasing):
        n = 10
        y = np.linspace(0, 1, n) if increasing else np.linspace(1, 0, n)
        x, w, r = isotonic_regression(y, increasing=increasing)
        assert_almost_equal(x, y)
        assert_almost_equal(r, np.arange(n + 1))

    def test_weights(self):
        w = np.array([1, 2, 5, 0.5, 0.5, 0.5, 1, 3])
        y = np.array([3, 2, 1, 10, 9, 8, 20, 10])
        x, wx, r = isotonic_regression(y, w)
        assert_almost_equal(x, [12/8, 12/8, 12/8, 9, 9, 9, 50/4, 50/4])
        assert_almost_equal(wx, [8, 1.5, 4])
        assert_almost_equal(r, [0, 3, 6, 8])

        # weights are like repeated observations, we repeat the 3rd element 5 times.
        w2 = np.array([1, 2, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 1, 3])
        y2 = np.array([3, 2, 1, 1, 1, 1, 1, 10, 9, 8, 20, 10])
        x2, wx2, r2 = isotonic_regression(y, w)
        assert_almost_equal(x2, x)
        assert_almost_equal(wx2, wx)
        assert_almost_equal(r2, r)

    def test_against_R_monotone(self):
        n = 100
        y = np.linspace(0, 1, num=n, endpoint=False)
        y = 5 * y + np.sin(10 * y)
        x, wx, r = isotonic_regression(y)

        # R code
        # library(monotone)
        # y <- 5 * ((1:n)-1)/n + sin(10 * ((1:n)-1)/n)
        # monotone(y)
        res = [
            0.0000000, 0.1498334, 0.2986693, 0.4455202, 0.5894183, 0.7294255,
            0.8646425, 0.9942177, 1.1173561, 1.2333269, 1.3414710, 1.4412074,
            1.5320391, 1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110,
            1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110,
            1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110,
            1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110,
            1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110,
            1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110, 1.5708110,
            1.5708110, 1.5708110, 1.5708110, 1.6241853, 1.7165453, 1.8177326,
            1.9272355, 2.0444597, 2.1687334, 2.2993145, 2.4353978, 2.5761233,
            2.7205845, 2.8678375, 3.0169106, 3.1668139, 3.3165492, 3.4651200,
            3.6115414, 3.7548499, 3.8941134, 4.0284398, 4.1569866, 4.2789690,
            4.3936679, 4.5004366, 4.5987081, 4.6880000, 4.7679197, 4.8381682,
            4.8656413, 4.8656413, 4.8656413, 4.8656413, 4.8656413, 4.8656413,
            4.8656413, 4.8656413, 4.8656413, 4.8656413, 4.8656413, 4.8656413,
            4.8656413, 4.8656413, 4.8656413, 4.8656413, 4.8656413, 4.8656413,
            4.8656413, 4.8656413, 4.8656413, 4.8656413,
        ]
        assert_almost_equal(x, res)

        # Test balance property: sum(y) == sum(x)
        assert_almost_equal(np.sum(x), np.sum(y))

        # Reverse order
        x, wx, r = isotonic_regression(-y, increasing=False)
        assert_almost_equal(-x, res)
