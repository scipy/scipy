import pytest
import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, assert_allclose
import scipy.special as sc


def test_ndtr():
    assert_equal(sc.ndtr(0), 0.5)
    assert_almost_equal(sc.ndtr(1), 0.84134474606)


class TestNdtri:

    def test_zero(self):
        assert sc.ndtri(0.5) == 0.0

    def test_asymptotes(self):
        assert_equal(sc.ndtri([0.0, 1.0]), [-np.inf, np.inf])

    def test_outside_of_domain(self):
        assert all(np.isnan(sc.ndtri([-1.5, 1.5])))


class TestNdtrDelta:

    @pytest.mark.parametrize('a, b', [(-50, 0), [0, 50]])
    def test_half(self, a, b):
        assert sc.ndtr_delta(a, b) == 0.5
        assert sc.ndtr_delta(b, a) == -0.5

    def test_zero(self):
        assert sc.ndtr_delta(3.0, 3.0) == 0

    # Expected values were computed with mpmath, e.g.:
    #
    #   >>> import mpmath
    #   >>> mpmath.mp.dps = 80
    #   >>> a = 10
    #   >>> b = 12
    #   >>> float(mpmath.ncdf(b) - mpmath.ncdf(a))
    #   7.619853022384043e-24
    #
    @pytest.mark.parametrize('a, b, expected',
                             [(10, 12, 7.619853022384043e-24),
                              (-1, 1.5, 0.7745375447996848),
                              (1, 5, 0.15865496727988518),
                              (-30, -25, 3.056696706382561e-138)])
    def test_values(self, a, b, expected):
        assert_allclose(sc.ndtr_delta(a, b), expected, rtol=1e-13)
        assert_allclose(sc.ndtr_delta(b, a), -expected, rtol=1e-13)


class TestTruncNdtr:

    def test_boundaries(self):
        assert sc.trunc_ndtr(1, 5, 1) == 0.0
        assert sc.trunc_ndtr(1, 5, 5) == 1.0

    def test_mid_symmetric(self):
        assert sc.trunc_ndtr(-3.2, 3.2, 0) == 0.5

    # Expected values were computed with mpmath, e.g.:
    #
    #   >>> import mpmath
    #   >>> mpmath.mp.dps = 80
    #   >>> a = 1
    #   >>> b = 5
    #   >>> x = 2
    #   >>> float((mpmath.ncdf(x) - mpmath.ncdf(a)) /
    #   ...       (mpmath.ncdf(b) - mpmath.ncdf(a)))
    #   0.8566080489842209
    #
    @pytest.mark.parametrize('x, a, b, expected',
                             [(2, 1, 5, 0.8566080489842209),
                              (1.000001, 1, 5, 1.5251372690210932e-06),
                              (11.25, 10, 15, 0.9999984803377245)])
    def test_values(self, x, a, b, expected):
        assert_allclose(sc.trunc_ndtr(a, b, x), expected, rtol=1e-10)
