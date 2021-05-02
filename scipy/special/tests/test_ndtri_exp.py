import math
import pytest
import numpy as np
from numpy.testing import assert_allclose

import scipy.special as sc


@pytest.fixture()
def uniform_random_points():
    def f(left, right, num_points=1000, seed=None):
        bit_generator = np.random.RandomState(seed)
        points = bit_generator.random_sample(num_points)
        return (right - left) * points + left
    return f


class TestNdtriExp(object):
    """Tests that ndtri_exp is sufficiently close to an inverse of log_ndtr

    Separate tests for the four intervals (-inf, -log(1e-6)),
    [-log(1e-6), -2), [-2, log(1 - exp(-2))), [log(1 - exp(-2)), inf).
    ndtri_exp(y) is computed in three different ways depending on if y 
    is in (-inf, -2), or one of the latter two intervals above. A separate
    test is provided for very small y.
    """
    @pytest.mark.parametrize('test_input',
                             [-15, -30, -1e2, -1e3, -1e6, -1e10, -1.7373e100,
                              -9.3493345e100, -5.8282223e300])
    def test_ndtri_exp_very_small_arg(self, test_input):
        y = test_input
        assert np.isclose(sc.log_ndtr(sc.ndtri_exp(y)), y, rtol=1e-12)

    def test_ndtri_exp_small_arg(self, uniform_random_points):
        points = uniform_random_points(1e-6, math.exp(-2),
                                       num_points=10**5, seed=561)
        assert_allclose(sc.ndtri_exp(sc.log_ndtr(points)), points,
                        rtol=1e-9)
                            
    def test_ndtri_exp_middle_arg(self, uniform_random_points):
        points = uniform_random_points(math.exp(-2), 1 - math.exp(-2),
                                       num_points=10**5, seed=1105)
        assert_allclose(sc.ndtri_exp(sc.log_ndtr(points)), points,
                        rtol=1e-14)

    def test_ndtri_exp_large_arg(self, uniform_random_points):
        points = uniform_random_points(1 - math.exp(-2), 1,
                                       num_points=10**5, seed=1729)
        assert_allclose(sc.ndtri_exp(sc.log_ndtr(points)), points,
                        rtol=1e-14)

    def test_ndtri_exp_outside_domain(self):
        assert np.isnan(sc.ndtri_exp(1.0))
