import numpy as np
import scipy.special as sc
import pytest
from numpy.testing import (assert_allclose, assert_array_equal,
                           suppress_warnings)


class TestNbdtr:
    def test(self):
        val = sc.nbdtr(0, 1, 0.5)
        assert_allclose(val, 0.5)

    def test_rounding(self):
        double_val = sc.nbdtr([0.1, 1.1, 2.1], 2, 0.5)
        int_val = sc.nbdtr([0, 1, 2], 2, 0.5)
        assert_array_equal(double_val, int_val)

    @pytest.mark.parametrize('k, n, p', [
        (np.inf, 2, 0.5),
        (1.0, np.inf, 0.5),
        (1.0, 2, np.inf)
    ])
    def test_inf(self, k, n, p):
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            val = sc.nbdtr(k, n, p)
        assert np.isnan(val)

    @pytest.mark.parametrize('k, n, p', [
        (-1.1, 1, 0.5),
        (2, 1, -0.1),
        (2, 1, 1.1),
    ])
    def test_domain(self, k, n, p):
        val = sc.nbdtr(k, n, p)
        assert np.isnan(val)


class TestNbdtrc:
    def test_value(self):
        val = sc.nbdtrc(0, 1, 0.5)
        assert_allclose(val, 0.5)

    def test_rounding(self):
        double_val = sc.nbdtrc([0.1, 1.1, 2.1], 2, 0.5)
        int_val = sc.nbdtrc([0, 1, 2], 2, 0.5)
        assert_array_equal(double_val, int_val)

    @pytest.mark.parametrize('k, n, p', [
        (np.inf, 2, 0.5),
        (1.0, np.inf, 0.5),
        (1.0, 2, np.inf)
    ])
    def test_inf(self, k, n, p):
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            val = sc.nbdtrc(k, n, p)
        assert np.isnan(val)

    @pytest.mark.parametrize('k, n, p', [
        (-1.1, 1, 0.5),
        (2, 1, -0.1),
        (2, 1, 1.1),
    ])
    def test_domain(self, k, n, p):
        val = sc.nbdtrc(k, n, p)
        assert np.isnan(val)

    def test_nbdtr_nbdtrc_sum_to_one(self):
        nbdtr_vals = sc.nbdtr([0, 1, 2], 2, 0.5)
        nbdtrc_vals = sc.nbdtrc([0, 1, 2], 2, 0.5)
        vals = nbdtr_vals + nbdtrc_vals
        assert_allclose(vals, [1.0, 1.0, 1.0])


class TestNbdtri:
    def test_value(self):
        val = sc.nbdtri(0, 1, 0.5)
        assert_allclose(val, 0.5)

    def test_rounding(self):
        double_val = sc.nbdtri([0.1, 1.1], 2, 0.5)
        int_val = sc.nbdtri([0, 1], 2, 0.5)
        assert_allclose(double_val, int_val)

    @pytest.mark.parametrize('k, n, p', [
        (np.inf, 2, 0.5),
        (1.0, np.inf, 0.5),
        (1.0, 2, np.inf)
    ])
    def test_inf(self, k, n, p):
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            val = sc.nbdtri(k, n, p)
        assert np.isnan(val)

    @pytest.mark.parametrize('k, n, p', [
        (-1.1, 1, 0.5),
        (2, 1, -0.1),
        (2, 1, 1.1),
    ])
    def test_domain(self, k, n, p):
        val = sc.nbdtri(k, n, p)
        assert np.isnan(val)

    def test_nbdtr_nbdtri_roundtrip(self):
        nbdtr_vals = sc.nbdtr([0, 1, 2], 2, 0.5)
        roundtrip_vals = sc.nbdtri([0, 1, 2], 2, nbdtr_vals)
        assert_allclose(roundtrip_vals, [0.5, 0.5, 0.5])
