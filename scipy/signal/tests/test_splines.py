# pylint: disable=missing-docstring
import numpy as np
from numpy import array
from numpy.testing import (assert_allclose, assert_array_equal,
                           assert_almost_equal, assert_raises)
import pytest
from pytest import raises

from scipy.signal._spline import (
    symiirorder1_ic, symiirorder2_ic_fwd, symiirorder2_ic_bwd)
import scipy.signal._splines as sp
from scipy import signal


class TestSymIIR:
    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64, np.complex64, np.complex128])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir1_ic(self, dtype, precision):
        c_precision = precision
        if precision == -1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Symmetrical initial conditions for a IIR filter of order 1 are:
        # x[0] + z1 * \sum{k = 0}^{n - 1} x[k] * z1^k

        # Check the initial condition for a low-pass filter
        # with coefficient b = 0.85 on a step signal. The initial condition is
        # a geometric series: 1 + b * \sum_{k = 0}^{n - 1} u[k] b^k.

        # Finding the initial condition corresponds to
        # 1. Computing the index n such that b**n < log(precision), which
        # corresponds to ceil(log(precision) / log(b))
        # 2. Computing the geometric series until n, this can be computed
        # using the partial sum formula: (1 - b**n) / (1 - b)
        # This holds due to the input being a step signal.
        b = 0.85
        n_exp = int(np.ceil(np.log(c_precision) / np.log(b)))
        expected = np.asarray([[(1 - b ** n_exp) / (1 - b)]], dtype=dtype)
        expected = 1 + b * expected

        # Create a step signal of size n + 1
        x = np.ones(n_exp + 1, dtype=dtype)
        assert_allclose(symiirorder1_ic(x, b, precision), expected,
                        atol=2e-6, rtol=2e-7)

        # Check the conditions for a exponential decreasing signal with base 2.
        # Same conditions hold, as the product of 0.5^n * 0.85^n is
        # still a geometric series
        b_d = array(b, dtype=dtype)
        expected = np.asarray(
            [[(1 - (0.5 * b_d) ** n_exp) / (1 - (0.5 * b_d))]], dtype=dtype)
        expected = 1 + b_d * expected

        # Create an exponential decreasing signal of size n + 1
        x = 2 ** -np.arange(n_exp + 1, dtype=dtype)
        assert_allclose(symiirorder1_ic(x, b, precision), expected,
                        atol=2e-6, rtol=2e-7)

    def test_symiir1_ic_fails(self):
        # Test that symiirorder1_ic fails whenever \sum_{n = 1}^{n} b^n > eps
        b = 0.85
        # Create a step signal of size 100
        x = np.ones(100, dtype=np.float64)

        # Compute the closed form for the geometrical series
        precision = 1 / (1 - b)
        assert_raises(ValueError, symiirorder1_ic, x, b, precision)

        # Test that symiirorder1_ic fails when |z1| >= 1
        assert_raises(ValueError, symiirorder1_ic, x, 1.0, -1)
        assert_raises(ValueError, symiirorder1_ic, x, 2.0, -1)
