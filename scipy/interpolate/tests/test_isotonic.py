import string

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
import pytest

from scipy.interpolate import IsotonicInterpolator


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

        assert_almost_equal(iso_str.y_, iso_num.y_)
        assert_almost_equal(iso_str(x_str), iso_num(x_num))

        if increasing:
            assert np.all(np.diff(iso_str(x_str)) >= 0)
        else:
            assert np.all(np.diff(iso_str(x_str)) <= 0)
