from math import sin

import pytest
from pytest import raises
from numpy.testing import assert_almost_equal

from scipy._lib._finite_differences import _derivative


class TestDerivative:
    def test_not_enough_points(self):
        with raises(ValueError, match="derivative order \\+ 1."):
            _derivative(lambda x: x, 0, order=5, n=3)

    def test_odd_n(self):
        with raises(ValueError, match="must be odd."):
            _derivative(lambda x: x, 0, n=2)

    @pytest.mark.parametrize(
        ["order", "n", "expected"],
        [(1, 3, 0.9800665778256222),
         (1, 5, 0.9800665778415817),
         (1, 7, 0.9800665778415817),
         (1, 9, 0.9800665778425575),
         (1, 11, 0.9800665778419693),
         (2, 3, -0.19866913669730477),
         (2, 5, -0.19866948364199996),
         (2, 7, -0.19866908899240918),
         (2, 9, -0.19866885805734644),
         (2, 11, -0.19866834834680008),
         (3, 5, -0.9992007221626407)]
    )
    def test_computation(self, order, n, expected):
        assert_almost_equal(_derivative(sin, 0.2, dx=1e-5, order=order, n=n),
                            expected)
