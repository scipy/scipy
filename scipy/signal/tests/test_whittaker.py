import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.signal import whittaker_henderson
from scipy.signal._whittaker import (
    _logdet_difference_matrix, _reml, _solveh_banded, _solve_WH_banded,
    _solve_WH_order2_fast,
)
from scipy.stats import special_ortho_group


def test_solveh_banded():
    n = 10
    # construct positive definite matrix A
    Q = special_ortho_group(n, seed=1234).rvs()
    A = Q @ np.diag(1 + np.arange(n)) @ Q.T
    b = np.arange(n)
    ab = np.zeros((n, n))
    ab[0, :] = np.diagonal(A, 0)
    for i in range(1, n):
        ab[i, :-i] = np.diagonal(A, i)
    x, logdet = _solveh_banded(ab, b, calc_logdet=True)
    assert_allclose(x, np.linalg.solve(A, b))
    assert_allclose(logdet, np.log(np.linalg.det(A)))

    # tridiagonal case
    x, logdet = _solveh_banded(ab[:2], b, calc_logdet=True)
    A_tri = (
        np.diag(np.diag(A)) + np.diag(np.diag(A, -1), -1) + np.diag(np.diag(A, 1), 1)
    )
    assert_allclose(logdet, np.log(np.linalg.det(A_tri)))


@pytest.mark.parametrize(
        ["order", "n"],
        [
            (1, 2), (1, 10),
            (2, 3), (2, 4), (2, 5), (2, 20),
            (5, 6), (5, 7), (5, 10), (5, 20),
            (6, 7), (6, 10), (6, 20),
        ])
def test_logdet_difference_matrix(order, n):
    p = order
    logdet = _logdet_difference_matrix(order=p, n=n)
    D = np.diff(np.eye(n), n=p, axis=0)  # shape (n-p, n)
    assert_allclose(logdet, np.log(np.linalg.det(D @ D.T)))
    eigvals = np.linalg.eigvals(D.T @ D)
    assert_allclose(logdet, np.sum(np.log(eigvals[eigvals > 1e-8])))


@pytest.mark.parametrize(
        ["signal", "lamb", "order", "weights", "msg"],
        [
            ([[1, 2, 3] * 3], 1, 2, None,
             "Input array signal must be of shape \\(n,\\)"),
            (np.zeros(2), 1, 2, None, "Input array signal must be at least of shape"),
            (np.arange(10), -0.9, 2, None, "Parameter lamb must be non-negative"),
            ([1, 2, 3], 1, 1.5, None,
             "Parameter order must be an integer larger equal 1."),
            ([1, 2, 3], 1, 0, None,
             "Parameter order must be an integer larger equal 1."),
            ([1, 2, 3], 1, 2, [0, 1],
             "Input array weights must have the same shape as the signal array."),
            ([1, 2, 3], 1, 2, [[0]],
             "Input array weights must have the same shape as the signal array."),
            ([1, 2, np.nan], 1, 2, None,
             "Input array signal must be finite."),
            ([1, 2, 3], 1, 2, [1, 2, np.nan],
             "Input array weights must be finite."),
            ([1, 2, np.nan], 1, 2, [1, 2, 3],
             "Input array weights must be zero for all non-finite"),
        ])
def test_whittaker_raises(signal, lamb, order, weights, msg):
    """Test that whittaker raises errors."""
    with pytest.raises(ValueError, match=msg):
        whittaker_henderson(signal, lamb=lamb, order=order, weights=weights)


def test_whittaker_small_data():
    """Test that whittaker works on a few data points."""
    # Should work on order + 1 data points. The first 2*order+1 are special.
    la = 1
    y = np.arange(2)
    x = whittaker_henderson(y, order=1, lamb=la)
    # Analytical solution for order=1 and n=2
    res = (np.array([[1 + la, la], [la, 1 + la]]) / (1 + 2 * la)) @ y
    assert_allclose(x, res, atol=1e-15)
    whittaker_henderson(np.arange(3), order=1, lamb=la)

    y = np.arange(3)
    x = whittaker_henderson(y, order=2, lamb=la)
    # Analytical solution for order=2 and n=3
    res = (
        np.array([
            [1 + 5 * la,     2 * la,        -la],
            [    2 * la, 1 + 2 * la,     2 * la],
            [       -la,     2 * la, 1 + 5 * la],
        ]) / (1 + 6 * la)) @ y
    assert_allclose(x, res, atol=1e-15)
    whittaker_henderson(np.arange(4), order=2, lamb=la)
    whittaker_henderson(np.arange(5), order=2, lamb=la)

    whittaker_henderson(np.arange(4), order=3, lamb=la)
    whittaker_henderson(np.arange(5), order=3, lamb=la)
    whittaker_henderson(np.arange(6), order=3, lamb=la)
    whittaker_henderson(np.arange(7), order=3, lamb=la)


@pytest.mark.parametrize("n", [3, 4, 5, 100])
def test_whittaker_direct_vs_fast_order2(n):
    """Test equivalent results"""
    rng = np.random.default_rng(42)
    signal = np.sin(2 * np.pi * np.linspace(0, 1, n)) + rng.standard_normal(n)
    x1, _ = _solve_WH_banded(signal, lamb=1.23, order=2)
    x2 = _solve_WH_order2_fast(signal, lamb=1.23)
    assert_allclose(x1, x2)


@pytest.mark.parametrize("order", [1, 2, 3])
def test_whittaker_limit_penalty(order):
    """Test that whittaker for close to zero and infinity penalty."""
    rng = np.random.default_rng(42)
    n = 100
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    noise = rng.standard_normal(n)
    signal = y + noise
    x = whittaker_henderson(signal, lamb=1e-7, order=order)
    assert_allclose(x, signal, rtol=1e-4, atol=1e-5)

    # In the limit of an infinite penalty, the smoothing results in an interpolation
    # polynom of degree = penalty order - 1 (order=2 => linear)
    x = whittaker_henderson(y, lamb=1e11, order=order)
    x_poly = np.arange(len(y))
    poly = np.polynomial.Polynomial.fit(x=x_poly, y=y, deg=order - 1)
    assert_allclose(x, poly(x_poly), rtol=10**(-6 + order), atol=1e-5)
    if order == 2:
        # Linear interpolation:
        # As the sine is positive fom 0 to pi and negative from pi to 2*pi, we expect
        # a negative slope.
        assert np.diff(x)[0] < 0


def test_whittaker_unpenalized():
    """Test whittaker for lamb=0."""
    n = 10
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    x = whittaker_henderson(y, lamb=0)
    assert_allclose(x, y)
    assert not np.may_share_memory(x, y)


def test_whittaker_weights():
    """Test that whittaker with weights of 1 is same as without weights."""
    rng = np.random.default_rng(42)
    n = 100
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    noise = rng.standard_normal(n)
    signal = y + noise
    x1 = whittaker_henderson(signal, lamb=1)
    w = np.ones_like(signal)
    x2 = whittaker_henderson(signal, lamb=1, weights=w)
    assert_allclose(x1, x2)

    # Multiplying penalty and weights by the same number does not change the result.
    x3 = whittaker_henderson(signal, lamb=3, weights=3*w)
    assert_allclose(x3, x1)


def test_whittaker_zero_weight_interpolation():
    """Test that whittaker interpolates where weights are zero."""
    n = 100
    signal = np.sin(2*np.pi * np.linspace(0, 1, n))
    signal[50:] += 2
    w = np.ones_like(signal)
    signal[40:60] = np.nan  # value does not matter, np.nan to check arg validation
    w[40:60] = 0
    # Note: interpolation is a polynomial of degree = 2 * order - 1.

    # order = 1 => linear interpolation
    x = whittaker_henderson(signal, lamb=1, order=1, weights=w)
    interp = np.interp(np.arange(40, 60), [39, 60], [x[39], x[60]])
    assert_allclose(x[40:60], interp)

    # order = 2 => cubic interpolation
    x = whittaker_henderson(signal, lamb=1, order=2, weights=w)
    poly = np.polynomial.Polynomial.fit(
        x=[38, 39, 60, 61], y=[x[38], x[39], x[60], x[61]], deg=3
    )
    assert_allclose(x[40:60], poly(np.arange(40, 60)))


@pytest.mark.parametrize("order", [1, 2, 3])
def test_reml_criterion(order):
    la = 1.2345
    n = 10
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    D = np.diff(np.eye(n), n=order, axis=0)
    M = D.T @ D

    def test_reml(la, y):
        A = np.eye(n) + la * M
        x = np.linalg.solve(A, y)
        resid = y - x
        return -0.5 * (resid @ resid + la * x @ M @ x
                       - np.log(np.linalg.det(la * D @ D.T))
                       + np.log(np.linalg.det(A)))

    r1 = _reml(lamb=la, y=y, order=order)
    r2 = test_reml(la=la, y=y)
    assert_allclose(r1, r2)
