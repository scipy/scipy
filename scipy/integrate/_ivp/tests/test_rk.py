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
    assert_allclose(np.sum(solver.B), 1, rtol=1e-15)
    assert_allclose(np.sum(solver.A, axis=1), solver.C, rtol=1e-14)


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
