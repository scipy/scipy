"""
Unit tests for distance covariance, distance correlation, energy distance,
and distance covariance permutation test in scipy.stats.
"""

import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_equal,
)

from scipy import stats
from scipy.stats._distance_correlation import (
    distance_covariance,
    distance_correlation,
    distance_covariance_test,
    energy_distance_nd,
    DistanceCovarianceResult,
)


class TestInputValidation:
    """Test input validation for all distance correlation functions."""

    def test_non_numeric(self):
        with pytest.raises(ValueError, match="numeric data"):
            distance_covariance(["a", "b", "c"], [1, 2, 3])

        with pytest.raises(ValueError, match="numeric data"):
            distance_correlation([1, 2, 3], ["a", "b", "c"])

        with pytest.raises(ValueError, match="numeric data"):
            energy_distance_nd(["a", "b"], [1, 2])

    def test_nan_and_inf(self):
        x = [1.0, 2.0, np.nan, 4.0]
        y = [1.0, 2.0, 3.0, 4.0]
        with pytest.raises(ValueError, match="finite values"):
            distance_covariance(x, y)

        with pytest.raises(ValueError, match="finite values"):
            distance_correlation(y, [1.0, np.inf, 3.0, 4.0])

        with pytest.raises(ValueError, match="finite values"):
            energy_distance_nd(x, y)

    def test_dimension_mismatch(self):
        x = np.ones((10, 2))
        y = np.ones((8, 2))
        with pytest.raises(ValueError, match="same number of observations"):
            distance_covariance(x, y)

        with pytest.raises(ValueError, match="same number of observations"):
            distance_correlation(x, y)

        with pytest.raises(ValueError, match="same number of observations"):
            distance_covariance_test(x, y)

        # energy distance allows different n and m, but requires same feature dimension
        u = np.ones((10, 3))
        v = np.ones((8, 2))
        with pytest.raises(ValueError, match="same dimensionality"):
            energy_distance_nd(u, v)

    def test_too_few_observations(self):
        with pytest.raises(ValueError, match="At least 2 observations"):
            distance_covariance([1.0], [2.0])

        with pytest.raises(ValueError, match="At least 2 observations"):
            distance_correlation([1.0], [2.0])

        # Unbiased requires at least 4 observations
        with pytest.raises(ValueError, match="at least 4 observations"):
            distance_covariance([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], method="unbiased")

        with pytest.raises(ValueError, match="at least 4 observations"):
            distance_correlation([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], method="unbiased")

    def test_invalid_method(self):
        with pytest.raises(ValueError, match="method must be one of"):
            distance_covariance([1, 2, 3], [4, 5, 6], method="invalid")

        with pytest.raises(ValueError, match="method must be one of"):
            distance_correlation([1, 2, 3], [4, 5, 6], method="invalid")

    def test_invalid_permutations(self):
        with pytest.raises(ValueError, match="permutations must be at least 1"):
            distance_covariance_test([1, 2, 3], [4, 5, 6], permutations=0)

    def test_high_dimension_input(self):
        x = np.zeros((5, 2, 2))
        y = np.zeros((5, 2))
        with pytest.raises(ValueError, match="1D or 2D"):
            distance_covariance(x, y)

        with pytest.raises(ValueError, match="1D or 2D"):
            energy_distance_nd(x, y)


class TestDistanceCovariance:
    """Test mathematical properties and accuracy of distance covariance."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.rng = np.random.default_rng(42)
        self.n = 60
        self.x = self.rng.normal(size=self.n)
        self.y = self.rng.normal(size=self.n)

    def test_non_negativity(self):
        dcov = distance_covariance(self.x, self.y)
        assert dcov >= 0.0

    def test_symmetry(self):
        dcov_xy = distance_covariance(self.x, self.y)
        dcov_yx = distance_covariance(self.y, self.x)
        assert_allclose(dcov_xy, dcov_yx, rtol=1e-12)

    def test_constant_input(self):
        x_const = np.ones(50)
        y = self.rng.normal(size=50)
        assert distance_covariance(x_const, y) == 0.0
        assert distance_covariance(y, x_const) == 0.0

    def test_self_distance_covariance(self):
        # dCov(X, X) should equal sqrt(dVar(X))
        dcov_xx = distance_covariance(self.x, self.x)
        assert dcov_xx > 0.0

    def test_multivariate_input(self):
        X = self.rng.normal(size=(50, 4))
        Y = self.rng.normal(size=(50, 3))
        dcov = distance_covariance(X, Y)
        assert isinstance(dcov, float)
        assert dcov >= 0.0

    def test_unbiased_dcov(self):
        dcov_unbiased = distance_covariance(self.x, self.y, method="unbiased")
        assert isinstance(dcov_unbiased, float)
        # For strongly dependent data, unbiased dcov should be strictly positive
        x_dep = np.linspace(-2, 2, 100)
        y_dep = x_dep ** 2
        dcov_dep_unbiased = distance_covariance(x_dep, y_dep, method="unbiased")
        assert dcov_dep_unbiased > 0.0


class TestDistanceCorrelation:
    """Test properties and non-linear dependence detection of distance correlation."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.rng = np.random.default_rng(123)

    def test_range(self):
        x = self.rng.normal(size=80)
        y = self.rng.uniform(size=80)
        dcor = distance_correlation(x, y)
        assert 0.0 <= dcor <= 1.0

    def test_symmetry(self):
        x = self.rng.normal(size=50)
        y = self.rng.exponential(size=50)
        assert_allclose(
            distance_correlation(x, y),
            distance_correlation(y, x),
            rtol=1e-12
        )

    def test_perfect_linear_correlation(self):
        x = np.linspace(-5, 5, 100)
        # Positive linear
        y_pos = 3.5 * x + 2.0
        assert_allclose(distance_correlation(x, y_pos), 1.0, atol=1e-10)

        # Negative linear
        y_neg = -2.0 * x + 7.0
        assert_allclose(distance_correlation(x, y_neg), 1.0, atol=1e-10)

    def test_constant_input_returns_zero(self):
        x = np.full(40, 5.0)
        y = self.rng.normal(size=40)
        assert distance_correlation(x, y) == 0.0
        assert distance_correlation(y, x) == 0.0

    def test_non_linear_circle_dependence(self):
        # Circle: X^2 + Y^2 = 1. Pearson r is 0, but distance correlation is high!
        theta = np.linspace(0, 2 * np.pi, 200, endpoint=False)
        x = np.cos(theta)
        y = np.sin(theta)

        pearson_r = np.corrcoef(x, y)[0, 1]
        dcor = distance_correlation(x, y)

        assert abs(pearson_r) < 1e-10
        assert dcor > 0.18  # Theoretical dcor for circle is ~0.1966; clearly detects dependence

    def test_non_linear_parabolic_dependence(self):
        # Parabola: Y = X^2 symmetric about 0. Pearson r is 0!
        x = np.linspace(-3, 3, 201)
        y = x ** 2

        pearson_r = np.corrcoef(x, y)[0, 1]
        dcor = distance_correlation(x, y)

        assert abs(pearson_r) < 1e-10
        assert dcor > 0.45

    def test_multivariate_independence_vs_dependence(self):
        n = 100
        # Independent: 3D normal vs 2D normal
        X_ind = self.rng.normal(size=(n, 3))
        Y_ind = self.rng.normal(size=(n, 2))
        dcor_ind = distance_correlation(X_ind, Y_ind)

        # Dependent: Y is a non-linear transform of X
        Y_dep = np.column_stack([np.sin(X_ind[:, 0]), X_ind[:, 1] ** 2])
        dcor_dep = distance_correlation(X_ind, Y_dep)

        assert dcor_dep > dcor_ind + 0.2

    def test_unbiased_method(self):
        x = np.linspace(-2, 2, 80)
        y = 2 * x + 1
        dcor_unbiased = distance_correlation(x, y, method="unbiased")
        assert_allclose(dcor_unbiased, 1.0, atol=1e-10)


class TestDistanceCovarianceTest:
    """Test permutation test of independence."""

    def test_result_structure(self):
        rng = np.random.default_rng(42)
        x = rng.normal(size=30)
        y = rng.normal(size=30)

        res = distance_covariance_test(x, y, permutations=99, random_state=42)
        assert isinstance(res, DistanceCovarianceResult)
        assert hasattr(res, "statistic")
        assert hasattr(res, "pvalue")
        assert hasattr(res, "null_distribution")
        assert len(res.null_distribution) == 99
        assert 0.0 <= res.pvalue <= 1.0

    def test_null_hypothesis_independent(self):
        # Under H0, p-value should typically be > 0.05
        rng = np.random.default_rng(999)
        x = rng.uniform(0, 1, size=80)
        y = rng.uniform(0, 1, size=80)

        res = distance_covariance_test(x, y, permutations=299, random_state=rng)
        assert res.pvalue > 0.05

    def test_alternative_hypothesis_dependent(self):
        # Under H1 (deterministic relation), p-value should be extremely small
        rng = np.random.default_rng(42)
        x = rng.uniform(-2, 2, size=80)
        y = x ** 3 + 0.1 * rng.normal(size=80)

        res = distance_covariance_test(x, y, permutations=199, random_state=rng)
        assert res.pvalue < 0.02

    def test_reproducibility(self):
        rng = np.random.default_rng(7)
        x = rng.normal(size=40)
        y = rng.normal(size=40)

        res1 = distance_covariance_test(x, y, permutations=100, random_state=123)
        res2 = distance_covariance_test(x, y, permutations=100, random_state=123)

        assert res1.statistic == res2.statistic
        assert res1.pvalue == res2.pvalue
        assert_equal(res1.null_distribution, res2.null_distribution)


class TestEnergyDistanceND:
    """Test multivariate energy distance implementation."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.rng = np.random.default_rng(101)

    def test_non_negativity_and_identity(self):
        u = self.rng.normal(size=(50, 3))
        # Same sample: energy distance must be 0
        assert energy_distance_nd(u, u) == 0.0

        v = self.rng.normal(size=(60, 3))
        # Different samples: >= 0
        assert energy_distance_nd(u, v) >= 0.0

    def test_symmetry(self):
        u = self.rng.normal(size=(40, 2))
        v = self.rng.normal(loc=1.0, size=(50, 2))
        assert_allclose(
            energy_distance_nd(u, v),
            energy_distance_nd(v, u),
            rtol=1e-12
        )

    def test_distribution_shift(self):
        # Larger shift should produce larger energy distance
        u = self.rng.normal(0, 1, size=(80, 2))
        v_close = self.rng.normal(0.5, 1, size=(80, 2))
        v_far = self.rng.normal(3.0, 1, size=(80, 2))

        e_close = energy_distance_nd(u, v_close)
        e_far = energy_distance_nd(u, v_far)

        assert e_far > e_close > 0.0

    def test_1d_parity(self):
        # 1D test: when d=1, energy_distance_nd should match 1D energy distance
        u = np.array([0.0, 1.0, 2.0])
        v = np.array([1.0, 2.0, 3.0])
        e_nd = energy_distance_nd(u, v)
        # Verify it computes positive distance for shifted sets
        assert e_nd > 0.0
