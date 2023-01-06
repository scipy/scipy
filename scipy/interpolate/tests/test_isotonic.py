import string

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
import pytest

from scipy.interpolate._pava_pybind import pava
from scipy.interpolate import isotonic_regression, IsotonicInterpolator


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

        # weights are like repeated observations, we repeat the 3rd element 5
        # times.
        w2 = np.array([1, 2, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 1, 3])
        y2 = np.array([3, 2, 1, 1, 1, 1, 1, 10, 9, 8, 20, 10])
        x2, wx2, r2 = isotonic_regression(y2, w2)
        assert_almost_equal(np.diff(x2[0:7]), 0)
        assert_almost_equal(x2[4:], x)
        assert_almost_equal(wx2, wx)
        assert_almost_equal(r2[1:] - 4, r[1:])

    def test_against_R_monotone(self):
        y = [0, 6, 8, 3, 5, 2, 1, 7, 9, 4]
        x, w, r = isotonic_regression(y)
        # R code
        # library(monotone)
        # options(digits=8)
        # monotone(c(0, 6, 8, 3, 5, 2, 1, 7, 9, 4))
        res = [
            0, 4.1666667, 4.1666667, 4.1666667, 4.1666667, 4.1666667,
            4.1666667, 6.6666667, 6.6666667, 6.6666667,
        ]
        assert_almost_equal(x, res)
        assert_equal(r, [0, 1, 7, 10])

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

        # Test increasing
        assert np.all(np.diff(x) >= 0)

        # Test balance property: sum(y) == sum(x)
        assert_almost_equal(np.sum(x), np.sum(y))

        # Reverse order
        x, wx, rinv = isotonic_regression(-y, increasing=False)
        assert_almost_equal(-x, res)
        assert_equal(rinv, r)


class TestIsotonicInterpolator:
    def test_permutation_invariance(self):
        # Check that the fitting (__init__) is permutation invariant.
        x = np.array([1, 2, 3, 4, 5, 6, 7])
        y = np.array([1, 41, 51, 1, 2, 5, 24])
        weights = np.array([1, 2, 3, 4, 5, 6, 7])
        iso1 = IsotonicInterpolator(x, y, weights=weights)

        rng = np.random.default_rng(314)
        perm = rng.permutation(x.shape[0])
        iso2 = IsotonicInterpolator(x[perm], y[perm], weights=weights[perm])

        assert_equal(iso1.x_, iso2.x_)
        assert_almost_equal(iso1.y_, iso2.y_)
        assert_almost_equal(iso1(x), iso2(x))

    @pytest.mark.parametrize("increasing", [True, False])
    def test_linspace(self, increasing):
        n = 10
        x = np.arange(n)
        y = np.linspace(0, 1, n) if increasing else np.linspace(1, 0, n)
        iso = IsotonicInterpolator(x, y, increasing=increasing)
        assert_almost_equal(iso(x), y)
        # Check linear interpolation
        assert_almost_equal(iso(x[:-1] + np.diff(x)/2), y[:-1] + np.diff(y)/2)
    
    @pytest.mark.parametrize("seed", list(range(310, 320)))
    @pytest.mark.parametrize("increasing", [True, False])
    def test_strings(self, seed, increasing):
        n = 10
        rng = np.random.default_rng(seed)
        x_str = np.array([string.ascii_lowercase[i] for i in range(n)])
        x_num = np.arange(n)
        y = np.arange(n)
        perm1, perm2 = rng.permutation(n), rng.permutation(n)
        idx_str = np.lexsort((y[perm2], x_str[perm1]))
        idx_num = np.lexsort((y[perm2], x_str[perm1]))
        assert_equal(idx_str, idx_num)

        iso_str = IsotonicInterpolator(
            x_str[perm1], y[perm2], increasing=increasing
        )
        iso_num = IsotonicInterpolator(
            x_num[perm1], y[perm2], increasing=increasing
        )

        assert not iso_str.x_is_numeric
        assert iso_num.x_is_numeric

        print(f"{iso_num.x_=}")
        print(f"{iso_str.x_=}")
        print(f"{iso_str.y_=}")
        assert_almost_equal(iso_str.y_, iso_num.y_)
        assert_almost_equal(iso_str(x_str), iso_num(x_num))
        
        if increasing:
            assert np.all(np.diff(iso_str(x_str)) >= 0)
        else:
            assert np.all(np.diff(iso_str(x_str)) <= 0)
