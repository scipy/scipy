import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal

import scipy.special as sc


def uniform_random_points(left, right, num_points=1000, seed=None):
    random_state = np.random.RandomState(seed)
    points = random_state.random_sample(num_points)
    return (right - left) * points + left


class TestNdtriExp(object):
    """Tests that ndtri_exp is sufficiently close to an inverse of log_ndtr

    Separate tests for the five intervals (-inf, -10),
    [-10, -2), [-2, -0.14542), [-0.14542, -1e-6), [-1e-6, 0).
    ndtri_exp(y) is computed in three different ways depending on if y 
    is in (-inf, -2), [-2, log(1 - exp(-2))], [log(1 - exp(-2), 0).
    Each of these intervals is given its own test with two additional tests
    for handling very small values and values very close to zero.
    """
    @pytest.mark.parametrize('test_input',
                             [-1e1, -1e2, -1e10, -1e20, -1e100])
    def test_very_small_arg(self, test_input):
        scale = test_input
        points = scale * uniform_random_points(0.5, 1, num_points=10**3,
                                               seed=561)
        assert_allclose(sc.log_ndtr(sc.ndtri_exp(points)), points,
                        rtol=1e-14)

    @pytest.mark.parametrize('test_input,expected_rtol',
                             [((-10, -2, 1105), 1e-14),
                              ((-2, -0.14542, 1729), 1e-12),
                              ((-0.14542, -1e-6, 2465), 1e-10),
                              ((-1e-6, 0, 2821), 1e-6)])
    def test_in_interval(self, test_input, expected_rtol):
        left, right, seed = test_input
        points = uniform_random_points(left, right, num_points=10**5,
                                       seed=seed)
        assert_allclose(sc.log_ndtr(sc.ndtri_exp(points)), points,
                        rtol=expected_rtol)

    def test_asymptotes(self):
        assert_equal(sc.ndtri_exp([-np.inf, -9e307, 0.0]),
                     [-np.inf, -np.inf, np.inf])

    def test_outside_domain(self):
        assert np.isnan(sc.ndtri_exp(1.0))


