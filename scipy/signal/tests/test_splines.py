# pylint: disable=missing-docstring
import numpy as np
from numpy import array
from numpy.testing import (assert_allclose, assert_array_equal,
                           assert_almost_equal, assert_raises)
import pytest
from pytest import raises

from scipy.signal._spline import (
    symiirorder1_ic, symiirorder2_ic_fwd, symiirorder2_ic_bwd)
from scipy.signal._splines import symiirorder1
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

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64, np.complex64, np.complex128])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir1(self, dtype, precision):
        c_precision = precision
        if precision == -1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Test for a low-pass filter with c0 = 0.15 and z1 = 0.85
        # using an unit step over 200 samples.
        c0 = 0.15
        z1 = 0.85
        n = 200
        signal = np.ones(n, dtype=dtype)

        # Find the initial condition. See test_symiir1_ic for a detailed
        # explanation
        n_exp = int(np.ceil(np.log(c_precision) / np.log(z1)))
        initial = np.asarray((1 - z1 ** n_exp) / (1 - z1), dtype=dtype)
        initial = 1 + z1 * initial

        # Forward pass
        # The transfer function for the system 1 / (1 - z1 * z^-1) when
        # applied to an unit step with initial conditions y0 is
        # 1 / (1 - z1 * z^-1) * (z^-1 / (1 - z^-1) + y0)

        # Solving the inverse Z-transform for the given expression yields:
        # y[n] = y0 * z1**n * u[n] +
        #        -z1 / (1 - z1) * z1**(k - 1) * u[k - 1] +
        #        1 / (1 - z1) * u[k - 1]
        # d is the Kronecker delta function, and u is the unit step

        # y0 * z1**n * u[n]
        pos = np.arange(n, dtype=dtype)
        comp1 = initial * z1**pos

        # -z1 / (1 - z1) * z1**(k - 1) * u[k - 1]
        comp2 = np.zeros(n, dtype=dtype)
        comp2[1:] = -z1 / (1 - z1) * z1**pos[:-1]

        # 1 / (1 - z1) * u[k - 1]
        comp3 = np.zeros(n, dtype=dtype)
        comp3[1:] = 1 / (1 - z1)

        expected_fwd = comp1 + comp2 + comp3

        # Reverse condition
        sym_cond = -c0 / (z1 - 1.0) * expected_fwd[-1]

        # Backward pass
        # The transfer function for the forward result is equivalent to
        # the forward system times c0 / (1 - z1 * z).

        # Computing a closed form for the complete expression is difficult
        # The result will be computed iteratively from the difference equation
        exp_out = np.zeros(n, dtype=dtype)
        exp_out[0] = sym_cond

        for i in range(1, n):
            exp_out[i] = c0 * expected_fwd[n - 1 - i] + z1 * exp_out[i - 1]

        exp_out = exp_out[::-1]

        out = symiirorder1(signal, c0, z1, precision)
        assert_allclose(out, exp_out, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir2_initial_fwd(self, dtype, precision):
        c_precision = precision
        if precision == -1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Compute the initial conditions for a order-two symmetrical low-pass
        # filter with r = 0.5 and omega = pi / 3 for an unit step input.
        r = np.asarray(0.5, dtype=dtype)
        omega = np.asarray(np.pi / 3.0, dtype=dtype)
        cs = 1 - 2 * r * np.cos(omega) + r**2

        # The index n for the initial condition is bound from 0 to the
        # first position where sin(omega * (n + 2)) = 0 => omega * (n + 2) = pi
        # For omega = pi / 3, the maximum initial condition occurs when
        # sqrt(3) / 2 * r**n < precision.
        # => n = log(2 * sqrt(3) / 3 * precision) / log(r)
        ub = np.ceil(np.log(c_precision / np.sin(omega)) / np.log(c_precision))
        lb = np.ceil(np.pi / omega) - 2
        n_exp = min(ub, lb)

        # The forward initial condition for a filter of order two is:
        # \frac{cs}{\sin(\omega)} \sum_{n = 0}^{N - 1} {
        #    r^(n + 1) \sin{\omega(n + 2)}} + cs
        # The closed expression for this sum is:
        # s[n] = 2 * r * np.cos(omega) -
        #        r**2 - r**(n + 2) * np.sin(omega * (n + 3)) / np.sin(omega) +
        #        r**(n + 3) * np.sin(omega * (n + 2)) / np.sin(omega) + cs
        fwd_initial_1 = (
            cs +
            2 * r * np.cos(omega) -
            r**2 -
            r**(n_exp + 2) * np.sin(omega * (n_exp + 3)) / np.sin(omega) +
            r**(n_exp + 3) * np.sin(omega * (n_exp + 2)) / np.sin(omega))

        # The second initial condition is given by
        # s[n] = 1 / np.sin(omega) * (
        #        r**2 * np.sin(3 * omega) -
        #        r**3 * np.sin(2 * omega) -
        #        r**(n + 3) * np.sin(omega * (n + 4)) +
        #        r**(n + 4) * np.sin(omega * (n + 3)))
        ub = np.ceil(np.log(c_precision / np.sin(omega)) / np.log(c_precision))
        lb = np.ceil(np.pi / omega) - 3
        n_exp = min(ub, lb)

        fwd_initial_2 = (
            cs + cs * 2 * r * np.cos(omega) +
            (r**2 * np.sin(3 * omega) -
             r**3 * np.sin(2 * omega) -
             r**(n_exp + 3) * np.sin(omega * (n_exp + 4)) +
             r**(n_exp + 4) * np.sin(omega * (n_exp + 3))) / np.sin(omega))

        expected = np.r_[fwd_initial_1, fwd_initial_2]

        n = 100
        signal = np.ones(n, dtype=dtype)

        out = symiirorder2_ic_fwd(signal, r, omega, precision)
        assert_allclose(out, expected, atol=4e-6, rtol=6e-7)
