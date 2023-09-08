import numpy as np
import pytest

from scipy.conftest import array_api_compatible
from scipy._lib._array_api import (
    _GLOBAL_CONFIG, array_namespace, as_xparray, copy, assert_equal
)


@pytest.mark.skipif(not _GLOBAL_CONFIG["SCIPY_ARRAY_API"],
        reason="Array API test; set environment variable SCIPY_ARRAY_API=1 to run it")
class TestArrayAPI:

    def test_array_namespace(self):
        x, y = np.array([0, 1, 2]), np.array([0, 1, 2])
        xp = array_namespace(x, y)
        assert 'array_api_compat.numpy' in xp.__name__

        _GLOBAL_CONFIG["SCIPY_ARRAY_API"] = False
        xp = array_namespace(x, y)
        assert 'array_api_compat.numpy' in xp.__name__
        _GLOBAL_CONFIG["SCIPY_ARRAY_API"] = True

    @array_api_compatible
    def test_asarray(self, xp):
        x, y = as_xparray([0, 1, 2], xp=xp), as_xparray(np.arange(3), xp=xp)
        ref = xp.asarray([0, 1, 2])
        assert_equal(x, ref)
        assert_equal(y, ref)

    @pytest.mark.filterwarnings("ignore: the matrix subclass")
    def test_raises(self):
        msg = "'numpy.ma.MaskedArray' are not supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace(np.ma.array(1), np.array(1))

        msg = "'numpy.matrix' are not supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace(np.array(1), np.matrix(1))

        msg = ("An argument was coerced to an object array, "
               "but object arrays are not supported.")
        with pytest.raises(TypeError, match=msg):
            array_namespace([object()])

    def test_array_likes(self):
        # should be no exceptions
        array_namespace([0, 1, 2])
        array_namespace(1, 2, 3)
        array_namespace(1)

    @array_api_compatible
    def test_copy(self, xp):
        for _xp in [xp, None]:
            x = xp.asarray([1, 2, 3])
            y = copy(x, xp=_xp)
            # with numpy we'd want to use np.shared_memory, but that's not specified
            # in the array-api
            x[0] = 10
            x[1] = 11
            x[2] = 12

            assert x[0] != y[0]
            assert x[1] != y[1]
            assert x[2] != y[2]
            assert id(x) != id(y)
