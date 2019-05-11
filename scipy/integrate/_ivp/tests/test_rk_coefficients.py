import pytest
from numpy.testing import assert_allclose
import numpy as np
from scipy.integrate import RK23, RK45
from scipy.integrate._ivp import dop853_coefficients


@pytest.mark.parametrize("solver", [RK23, RK45])
def test_coefficient_properties(solver):
    assert_allclose(np.sum(solver.B), 1, rtol=1e-15)
    assert_allclose(np.sum(solver.A, axis=1), solver.C, rtol=1e-15)


def test_coefficient_properties_dop853():
    assert_allclose(np.sum(dop853_coefficients.B), 1, rtol=1e-15)
    assert_allclose(np.sum(dop853_coefficients.A, axis=1),
                    dop853_coefficients.C,
                    rtol=1e-14)
