import numpy as np
from numpy.testing import assert_allclose
import pytest

from scipy._lib._bounds import _validate_bounds


def test_bounds():
    l_bounds = np.array([5., 0.])
    u_bounds = np.array([5., 0.95])
    x0 = np.array([5., 0.52])

    bounds = _validate_bounds(l_bounds=l_bounds, u_bounds=u_bounds, x0=x0)
    assert_allclose(bounds[0], l_bounds)
    assert_allclose(bounds[1], u_bounds)

    u_bounds = np.array([5.])
    bounds = _validate_bounds(l_bounds, u_bounds, x0=x0)
    assert_allclose(bounds[0], l_bounds)
    assert_allclose(bounds[1], np.array([5., 5.]))


def test_raises():
    sample = np.array([[0, 0], [1, 1], [0.5, 0.5]])

    msg = "An upper bound is less than the corresponding lower bound"
    with pytest.raises(ValueError, match=msg):
        bounds = np.array([[-2, 6], [6, 5]])
        _validate_bounds(l_bounds=bounds[0], u_bounds=bounds[1], x0=sample)

    msg = "The number of bounds is not compatible"
    with pytest.raises(ValueError, match=msg):
        l_bounds, u_bounds = [-2, 0, 2], [6, 5]
        _validate_bounds(l_bounds=l_bounds, u_bounds=u_bounds, x0=sample)

    with pytest.raises(ValueError, match=msg):
        bounds = np.array([[-2, 0, 2], [6, 5, 5]])
        _validate_bounds(l_bounds=bounds[0], u_bounds=bounds[1], x0=sample)
