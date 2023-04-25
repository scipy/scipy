import numpy as np
from numpy.testing import assert_equal
import pytest

from scipy.conftest import array_api_compatible
from scipy._lib._array_api import (
    USE_ARRAY_API, namespace_from_arrays, asarray, asarray_namespace
)


if not USE_ARRAY_API:
    pytest.skip(
        "Array API test; set environment variable array_api_dispatch=1 to run it",
        allow_module_level=True
    )


def test_namespace_from_arrays():
    x, y = [0, 1, 2], np.arange(3)
    xp = namespace_from_arrays(x, y)
    assert xp.__name__ == 'array_api_compat.numpy'


@array_api_compatible
def test_asarray(xp):
    x, y = asarray([0, 1, 2], xp=xp), asarray(np.arange(3), xp=xp)
    ref = np.array([0, 1, 2])
    assert_equal(x, ref)
    assert_equal(y, ref)


@array_api_compatible
def test_asarray_namespace(xp):
    x, y = [0, 1, 2], xp.arange(3)
    x, y, xp_ = asarray_namespace(x, y)
    assert xp_.__name__ == 'array_api_compat.numpy'
    ref = np.array([0, 1, 2])
    assert_equal(x, ref)
    assert_equal(y, ref)
    assert type(x) == type(y)
