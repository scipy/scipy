import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from scipy.special import pow1p


class TestPow1p:

    @pytest.mark.parametrize('x, y', [
        (1.2, 3.4),        # "ordinary" value
        (-2.0, 3.0),       # (1+x) < 0
        (1e-37, 1e38),     # normal float32, but < 2e-16
        (-1e-45, 3.4e38),  # subnormal float32, and < 2e-16
    ])
    def test_float32(self, x, y):
        # float32 inputs should give the same answer as float64 inputs
        x32, y32 = np.float32(x), np.float32(y)
        x64, y64 = np.float64(x32), np.float64(y32)
        res32 = pow1p(x32, y32)
        assert isinstance(res32, np.float32)
        res64 = pow1p(x64, y64)
        assert isinstance(res64, np.float64)
        assert_equal(res32, np.float32(res64))

    @pytest.mark.parametrize('rule, x, y, expected', [
        (
            "pow(anything, +0) is 1",
            [1.5, -1.0, -np.inf, np.inf, np.nan], 0.0, 1.0
        ),
        (
            "pow(anything, -0) is 1",
            [1.5, -1.0, -np.inf, np.inf, np.nan], -0.0, 1.0
        ),
        (
            "pow(±0, negative odd integer) is ±∞ and signals divideByZero",
            -1, [-3.0, -1.0], np.inf
        ),
        (
            "pow(±0, -∞) is +∞ with no exception",
            -1.0, -np.inf, np.inf
        ),
        (
            "pow(±0, +∞) is +0 with no exception",
            -1.0, np.inf, 0.0
        ),
        (
            "pow(±0, positive odd integer) is ±0",
            -1.0, [1.0, 3.0], 0.0,
        ),
        (
            "pow(-1, ±∞) is 1 with no exception",
            -2.0, [np.inf, -np.inf], 1.0
        ),
        (
            "pow(+1, anything) is 1",
            0.0, [1.5, -1.0, -np.inf, np.inf, np.nan], 1.0
        ),
        (
            "pow(|x| < 1, +∞) is +0",
            [np.nextafter(-2, 0), -1.0, np.nextafter(0, -1)], np.inf, 0.0
        ),
        (
            "pow(|x| > 1, +∞) is +∞",
            [np.nextafter(-2, -3), np.nextafter(0, 1), -np.inf, np.inf],
            np.inf, np.inf
        ),
        (
            "pow(|x| < 1, -∞) is +∞",
            [np.nextafter(-2, 0), -1.0, np.nextafter(0, -1)], -np.inf, np.inf
        ),
        (
            "pow(|x| > 1, -∞) is +0",
            [np.nextafter(-2, -3), np.nextafter(0, 1), -np.inf, np.inf],
            -np.inf, 0.0
        ),
        (
            "pow(+∞, y < 0) is +0",
            np.inf, [-np.inf, -2.5, np.nextafter(0, -1)], 0.0
        ),
        (
            "pow(+∞, y > 0) is +∞",
            np.inf, [np.inf, +2.5, np.nextafter(0, 1)], np.inf
        ),
        (
            "pow(-∞, negative odd integer y) is -0",
            -np.inf, [-3.0, -1.0], -0.0
        ),
        (
            "pow(-∞, positive odd integer y) is -∞",
            -np.inf, [1.0, 3.0], -np.inf
        ),
        (
            "pow(-∞, negative finite y not an odd integer) is +0",
            -np.inf, [-2.0, -3.5], +0.0,
        ),
        (
            "pow(-∞, positive finite y not an odd integer) is +∞",
            -np.inf, [2.0, 3.5], np.inf,
        ),
        (
            "pow(±0, negative finite y not an odd integer) is +∞ "
            "and signals divideByZero",
            -1, [-2.0, -3.5], np.inf,
        ),
        (
            "pow(±0, positive finite y not an odd integer) is +0",
            -1, [2.0, 3.5], 0,
        ),
        (
            "pow(finite negative x, finite non-integer y) "
            "signals the invalid operation exception",
            -2.5, [1.5, -2.5], np.nan,
        )
    ])
    def test_ieee754_special_values(self, rule, x, y, expected):
        # Special value inputs should conform to IEEE-754 spec of pow().
        assert_equal(pow1p(x, y), expected, err_msg=rule)

    @pytest.mark.parametrize('rule, x, y, expected', [
        ("spurious underflow", -1e-16, 7e18, 9.859676543759570e-305),
        ("spurious underflow & overflow 1", -1e-16, 1e20, 0.0),
        ("spurious underflow & overflow 2", 1e-16, 1e20, np.inf),
    ])
    def test_decomposition(self, rule, x, y, expected):
        # Test that the decomposition 1+x == s+t is done properly to avoid
        # spurious underflow and overflow.  This amounts to ensuring that
        # t has the same sign as x.  The 'expected' value is computed using
        # mpmath.
        if expected == 0 or expected == np.inf:
            assert_equal(pow1p(x, y), expected, err_msg=rule)
        else:
            assert_allclose(pow1p(x, y), expected, rtol=1e-15, err_msg=rule)
