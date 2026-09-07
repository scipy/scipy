import numpy as np
import pytest
import scipy.special as sc


class TestRgamma:

    def test_gh_11315(self):
        assert sc.rgamma(-35) == 0

    def test_rgamma_zeros(self):
        x = np.array([0, -10, -100, -1000, -10000])
        assert np.all(sc.rgamma(x) == 0)


class TestLogGammainc:
    @pytest.mark.parametrize("a,z", [
        (0, 1),      # a must be > 0
        (-1, 1),     # a negative
        (1, -1),     # z negative
        (0, 0),      # a = 0
    ])
    def test_domain_error_returns_nan(self, a, z):
        assert np.isnan(sc.log_gammainc(a, z))

    def test_both_inf_returns_nan(self):
        assert np.isnan(sc.log_gammainc(np.inf, np.inf))

    def test_z_zero_returns_neg_inf(self):
        assert sc.log_gammainc(1.0, 0.0) == -np.inf

    @pytest.mark.parametrize("a", [0.5, 1.0, 2.0, 100.0])
    def test_z_zero_returns_neg_inf_parametrized(self, a):
        assert sc.log_gammainc(a, 0.0) == -np.inf

    def test_z_inf_returns_zero(self):
        assert sc.log_gammainc(1.0, np.inf) == 0.0

    @pytest.mark.parametrize("a", [0.5, 1.0, 2.0, 100.0])
    def test_z_inf_returns_zero_parametrized(self, a):
        assert sc.log_gammainc(a, np.inf) == 0.0

    def test_a_inf_finite_z_returns_neg_inf(self):
        assert sc.log_gammainc(np.inf, 1.0) == -np.inf

    @pytest.mark.parametrize("z", [0.5, 1.0, 2.0, 100.0])
    def test_a_inf_finite_z_returns_neg_inf_parametrized(self, z):
        assert sc.log_gammainc(np.inf, z) == -np.inf


class TestLogGammaincc:
    """Tests for edge cases of log_gammaincc(a, z) = log(Q(a, z)) = log(1 - P(a, z))."""

    @pytest.mark.parametrize("a,z", [
        (0, 1),
        (-1, 1),
        (1, -1),
        (0, 0),
    ])
    def test_domain_error_returns_nan(self, a, z):
        assert np.isnan(sc.log_gammaincc(a, z))

    def test_both_inf_returns_nan(self):
        assert np.isnan(sc.log_gammaincc(np.inf, np.inf))

    def test_z_zero_returns_zero(self):
        assert sc.log_gammaincc(1.0, 0.0) == 0.0

    @pytest.mark.parametrize("a", [0.5, 1.0, 2.0, 100.0])
    def test_z_zero_returns_zero_parametrized(self, a):
        assert sc.log_gammaincc(a, 0.0) == 0.0

    def test_z_inf_returns_neg_inf(self):
        assert sc.log_gammaincc(1.0, np.inf) == -np.inf

    @pytest.mark.parametrize("a", [0.5, 1.0, 2.0, 100.0])
    def test_z_inf_returns_neg_inf_parametrized(self, a):
        assert sc.log_gammaincc(a, np.inf) == -np.inf

    def test_a_inf_finite_z_returns_zero(self):
        assert sc.log_gammaincc(np.inf, 1.0) == 0.0

    @pytest.mark.parametrize("z", [0.5, 1.0, 2.0, 100.0])
    def test_a_inf_finite_z_returns_zero_parametrized(self, z):
        assert sc.log_gammaincc(np.inf, z) == 0.0
