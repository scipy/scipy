import numpy as np
import pytest
from scipy.interpolate import BSpline, _dierckx
from scipy.interpolate._bsplines import _periodic_knots
from scipy.interpolate._fitpack_repro import Fperiodic


@pytest.mark.parametrize("k", [0, 2, 5])
@pytest.mark.parametrize("mode", ["regular", "repeated_skip", "knots_and_jumps"])
@pytest.mark.parametrize("dim", [1, 3])
@pytest.mark.parametrize("weight_mode", ["unit", "nonuniform", "zero"])
def test_fpback_residual_matches_spline_evaluation(k, mode, dim, weight_mode):
    rng = np.random.default_rng(4687 + k + dim)
    t = np.r_[np.repeat(0., k + 1), np.arange(1., 10.),
              np.repeat(10., k + 1)]
    nc = len(t) - k - 1
    patterns = {
        "regular": np.linspace(0., 10., 30),
        "repeated_skip": np.repeat(
            np.array([0., .2, 2.5, 2.5, 3.8, 9., 10.]), [4, 4, 4, 4, 4, 4, 6]
        ),
        "knots_and_jumps": np.repeat(np.array([0., 1., 1., 5., 9., 10.]), 5),
    }
    x = patterns[mode]
    R = np.zeros((len(x), k + 1))
    R[:nc, 0] = 1.
    coefficients = rng.normal(size=(nc, dim))
    yw = np.zeros((len(x), dim))
    yw[:nc] = coefficients
    y = np.ascontiguousarray(rng.normal(size=(len(x), dim)))
    w = np.ones(len(x)) if weight_mode == "unit" else rng.uniform(.1, 2., len(x))
    if weight_mode == "zero":
        w[::7] = 0.

    c, residuals, fp = _dierckx.fpback(R, nc, x, y, t, k, w, yw)
    expected = np.sum((w[:, None] * (BSpline(t, c, k)(x) - y)) ** 2, axis=1)
    np.testing.assert_allclose(c, coefficients, rtol=0, atol=0)
    np.testing.assert_allclose(residuals, expected, rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(fp, expected.sum(), rtol=1e-11, atol=1e-11)


@pytest.mark.parametrize("p", [0.1, 1., 10.])
def test_periodic_smoothing_residual_matches_spline(p):
    # Repeated observations exercise knot-interval tracking.
    k = 3
    x = np.array([0., .1, .2, 2.5, 2.5, 3.8, 4.1, 5., 6., 7., 8., 9., 10.])
    y = np.sin(x)[:, None]
    w = np.linspace(.5, 1.5, len(x))
    t = _periodic_knots(np.linspace(0., 10., 11), k)
    objective = Fperiodic(x, y, t, k=k, s=.2, w=w)
    actual = objective(p)
    values = objective.spl(x[:-1])
    expected = np.sum((w[:-1, None] * (values - y[:-1])) ** 2) - .2
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("mode", ["unique", "repeated"])
@pytest.mark.parametrize("s", [.05, .2, 1.])
def test_periodic_splrep_achieves_requested_residual(mode, s):
    from scipy.interpolate import make_splrep

    x = np.array([0., .1, .2, 2.5, 2.6, 3.8, 4.1, 5., 6., 7., 8., 9., 10.])
    if mode == "repeated":
        x[4] = x[3]
    y = np.sin(x) + .2 * np.cos(2 * x)
    w = np.linspace(.5, 1.5, len(x))
    t = _periodic_knots(np.linspace(0., 10., 11), 3)
    spline = make_splrep(x, y, w=w, k=3, s=s, t=t, bc_type="periodic")

    # Periodic fitting omits the last sample from the residual.
    fp = np.sum((w[:-1] * (spline(x[:-1]) - y[:-1])) ** 2)
    np.testing.assert_allclose(fp, s, rtol=1e-3, atol=0)


@pytest.mark.parametrize("k", [0, 2, 5])
def test_fpback_residual_with_extrapolation(k):
    t = np.r_[np.repeat(0., k + 1), np.arange(1., 10.),
              np.repeat(10., k + 1)]
    nc = len(t) - k - 1
    x = np.repeat(np.array([-1., 0., .2, 2.5, 2.5, 7., 11.]),
                  [4, 4, 4, 4, 4, 4, 6])
    R = np.zeros((len(x), k + 1))
    R[:nc, 0] = 1.
    yw = np.zeros((len(x), 1))
    yw[:nc, 0] = np.linspace(-1., 1., nc)
    y = np.sin(x)[:, None]
    w = np.linspace(.5, 1.5, len(x))

    c, residuals, fp = _dierckx.fpback(R, nc, x, y, t, k, w, yw, True)
    expected = (w * (BSpline(t, c[:, 0], k)(x) - y[:, 0])) ** 2
    np.testing.assert_allclose(residuals, expected, rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(fp, expected.sum(), rtol=1e-11, atol=1e-11)
