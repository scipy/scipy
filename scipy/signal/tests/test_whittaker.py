import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.signal import whittaker_henderson
from scipy.signal._whittaker import _solve_WH_banded, _solve_WH_order2_fast


@pytest.mark.parametrize(
        ["signal", "lamb", "msg"],
        [
            ([[1, 2, 3] * 3], 1, "Input signal array must be of shape \\(n,\\)"),
            (np.zeros(2), 1, "Input signal array must be at least of shape"),
            (np.arange(10), -0.9, "Parameter lamb must be non-negative"),
        ])
def test_whittaker_raises(signal, lamb, msg):
    """Test that whittaker raises errors."""
    with pytest.raises(ValueError, match=msg):
        whittaker_henderson(signal, lamb=lamb)


def test_whittaker_small_data():
    """Test that whittaker works on a few data points."""
    # Should work on order + 1 data points (so far, order=2 always)
    whittaker_henderson(np.zeros(3))
    whittaker_henderson(np.zeros(4))
    whittaker_henderson(np.zeros(5))


@pytest.mark.parametrize("n", [3, 4, 5, 100])
def test_whittaker_direct_vs_fast_order2(n):
    """Test equivalent results"""
    rng = np.random.default_rng(42)
    signal = np.sin(2 * np.pi * np.linspace(0, 1, n)) + rng.standard_normal(n)
    x1 = _solve_WH_banded(signal, lamb=1.23)
    x2 = _solve_WH_order2_fast(signal, lamb=1.23)
    assert_allclose(x1, x2)


def test_whittaker_order_2():
    """Test that whittaker reproduces a signal."""
    rng = np.random.default_rng(42)
    n = 100
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    noise = rng.standard_normal(n)
    signal = y + noise
    x = whittaker_henderson(signal, lamb=1e-6)
    assert_allclose(x, signal, rtol=1e-4, atol=1e-5)

    # In the limit of an infinite penalty, the smoothing results in a linear
    # interpolation.
    x = whittaker_henderson(y, lamb=1e9)
    assert_allclose(np.diff(x, n=2), 0, atol=1e-6)
    # As the sine is positive fom 0 to pi and negative from pi to 2*pi, we expect
    # a negative slope.
    assert np.diff(x)[0] < 0
