import numpy as np
import pytest

from scipy._lib._bounds import _validate_bounds


def test_bounds():
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
