import numpy as np
from numpy.testing import assert_equal
import pytest

from scipy.conftest import array_api_compatible
from scipy._lib._array_api import (
    _GLOBAL_CONFIG, array_namespace, as_xparray,
)


if not _GLOBAL_CONFIG["SCIPY_ARRAY_API"]:
    pytest.skip(
        "Array API test; set environment variable SCIPY_ARRAY_API=1 to run it",
        allow_module_level=True
    )


def to_numpy(array, xp):
    """Convert `array` into a NumPy ndarray on the CPU. From sklearn."""
    xp_name = xp.__name__

    if xp_name in {"array_api_compat.torch", "torch"}:
        return array.cpu().numpy()
    elif xp_name == "cupy.array_api":
        return array._array.get()
    elif xp_name in {"array_api_compat.cupy", "cupy"}:  # pragma: nocover
        return array.get()

    return np.asarray(array)


def test_array_namespace():
    x, y = np.array([0, 1, 2]), np.array([0, 1, 2])
    xp = array_namespace(x, y)
    assert 'array_api_compat.numpy' in xp.__name__

    _GLOBAL_CONFIG["SCIPY_ARRAY_API"] = False
    xp = array_namespace(x, y)
    assert 'array_api_compat.numpy' in xp.__name__
    _GLOBAL_CONFIG["SCIPY_ARRAY_API"] = True


@array_api_compatible
def test_asarray(xp):
    x, y = as_xparray([0, 1, 2], xp=xp), as_xparray(np.arange(3), xp=xp)
    ref = np.array([0, 1, 2])
    assert_equal(x, ref)
    assert_equal(y, ref)


@array_api_compatible
def test_to_numpy(xp):
    x = xp.asarray([0, 1, 2])
    x = to_numpy(x, xp=xp)
    assert isinstance(x, np.ndarray)


@pytest.mark.filterwarnings("ignore: the matrix subclass")
def test_raises():
    msg = "'numpy.ma.MaskedArray' are not supported"
    with pytest.raises(TypeError, match=msg):
        array_namespace(np.ma.array(1), np.array(1))

    msg = "'numpy.matrix' are not supported"
    with pytest.raises(TypeError, match=msg):
        array_namespace(np.array(1), np.matrix(1))

    msg = "Only support Array API"
    with pytest.raises(TypeError, match=msg):
        array_namespace([0, 1, 2])

    with pytest.raises(TypeError, match=msg):
        array_namespace(1)
