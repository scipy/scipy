import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.signal import whittaker_handerson
from scipy.signal._whittaker import _solve_WH_order_direct, _solve_WH_order2_fast


def test_whittaker_small_data():
    """Test that whittaker works on a few data points."""
    msg = "Input signal array must be at least of shape*"
    with pytest.raises(ValueError, match=msg):
        whittaker_handerson(np.zeros(2))
    # Should work on order + 1 data points (so far, order=2 always)
    whittaker_handerson(np.zeros(3))
    whittaker_handerson(np.zeros(4))
    whittaker_handerson(np.zeros(5))


def test_whittaker_direct_vs_fast_order2():
    """Test equivalent results"""
    rng = np.random.default_rng(42)
    n = 100
    signal = np.sin(2 * np.pi * np.linspace(0, 1, n)) + rng.standard_normal(n)
    x1 = _solve_WH_order_direct(signal, lamb=1/n)
    x2 = _solve_WH_order2_fast(signal, lamb=1/n)
    assert_allclose(x1, x2)


def test_whittaker_order_2():
    """Test that whittaker reproduces a signal."""
    rng = np.random.default_rng(42)
    n = 100
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    noise = rng.standard_normal(n)
    signal = y + noise
    x = whittaker_handerson(signal, lamb=1e-6)
    assert_allclose(x, signal, rtol=1e-4, atol=1e-5)

    # In the limit of an infinite penalty, the smoothing results in a linear
    # interpolation.
    x = whittaker_handerson(y, lamb=1e9)
    assert_allclose(np.diff(x, n=2), 0, atol=1e-6)
    # As the sine is positive fom 0 to pi and negative from pi to 2*pi, we expect
    # a negative slope.
    assert np.diff(x)[0] < 0
