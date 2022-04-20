import pytest
from numpy.testing import assert_allclose, assert_
import numpy as np
from scipy.integrate import (
    RK23, RK45, DOP853, Tsit5, Verner65, Verner76, Verner87, Verner98)
from scipy.integrate._ivp import coefs

rk_methods = [
    RK23, RK45, DOP853, Tsit5, Verner65, Verner76, Verner87, Verner98
]
coef_sources = [
    RK23, RK45, DOP853, coefs.tsit5,
    coefs.verner65, coefs.verner76, coefs.verner87, coefs.verner98
]


@pytest.mark.parametrize("solver", coef_sources)
def test_coefficient_properties(solver):
    # the coefficients can only be expect to satisfy these properties to the
    # extent that they are of the same order of magnitude, which determines
    # the roundoff error (i.e., catastrophic cancellation)

    absB = np.abs(solver.B)
    rtol_B = np.max(absB) / np.min(absB[absB > 0]) * 6e-16
    assert_allclose(np.sum(solver.B), 1, rtol=max(1e-15, rtol_B))

    for i, a in enumerate(solver.A[1:], start=1):
        absA = np.abs(a[:i])
        rtol_A = np.max(absA) / np.min(absA[absA > 0]) * 6e-16
        print(i, solver.__name__)
        assert_allclose(np.sum(a), solver.C[i], rtol=max(1e-14, rtol_A))


@pytest.mark.parametrize("solver_class", rk_methods)
def test_error_estimation(solver_class):
    step = 0.2
    solver = solver_class(lambda t, y: y, 0, [1], 1, first_step=step)
    solver.step()
    error_estimate = solver._estimate_error(solver.K, step)
    error = solver.y - np.exp([step])
    assert_(np.abs(error) < np.abs(error_estimate))


@pytest.mark.parametrize("solver_class", rk_methods)
def test_error_estimation_complex(solver_class):
    h = 0.2
    solver = solver_class(lambda t, y: 1j * y, 0, [1j], 1, first_step=h)
    solver.step()
    err_norm = solver._estimate_error_norm(solver.K, h, scale=[1])
    assert np.isrealobj(err_norm)
