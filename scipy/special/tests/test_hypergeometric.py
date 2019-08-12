import pytest
import numpy as np
from numpy.testing import assert_allclose
from numpy.testing import assert_equal

import scipy.special as sc


class TestHyperu:

    def test_negative_x(self):
        a, b, x = np.meshgrid(
            [-1, -0.5, 0, 0.5, 1],
            [-1, -0.5, 0, 0.5, 1],
            np.linspace(-100, -1, 10),
        )
        assert np.all(np.isnan(sc.hyperu(a, b, x)))

    def test_special_cases(self):
        assert sc.hyperu(0, 1, 1) == 1.0


class TestHyp1f1:

    @pytest.mark.parametrize('a, b, x', [
        (np.nan, 1, 1),
        (1, np.nan, 1),
        (1, 1, np.nan)
    ])
    def test_nan_inputs(self, a, b, x):
        assert np.isnan(sc.hyp1f1(a, b, x))

    def test_poles(self):
        assert_equal(sc.hyp1f1(1, [0, -1, -2, -3, -4], 0.5), np.infty)

    @pytest.mark.parametrize('a, b, x, result', [
        (-1, 1, 0.5, 0.5),
        (1, 1, 0.5, 1.6487212707001281468),
        (2, 1, 0.5, 2.4730819060501922203),
        (1, 2, 0.5, 1.2974425414002562937),
        (-10, 1, 0.5, -0.38937441413785204475)
    ])
    def test_special_cases(self, a, b, x, result):
        # Hit all the special case branches at the beginning of the
        # function. Desired answers computed using Mpmath.
        assert_allclose(sc.hyp1f1(a, b, x), result, atol=0, rtol=1e-15)

    def test_gh_3492(self):
        desired = 0.99973683897677527773  # Computed using Mpmath
        assert_allclose(
            sc.hyp1f1(0.01, 150, -4),
            desired,
            atol=0,
            rtol=1e-15
        )

    def test_gh_3593(self):
        desired = 1.0020033381011970966  # Computed using Mpmath
        assert_allclose(
            sc.hyp1f1(1, 5, 0.01),
            desired,
            atol=0,
            rtol=1e-15
        )
