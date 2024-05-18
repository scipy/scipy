import numpy as np
import pytest

from scipy.conftest import array_api_compatible
from scipy._lib._array_api import (
    _GLOBAL_CONFIG, array_namespace, _asarray, copy, xp_assert_equal, is_numpy,
    xp_mean, xp_assert_close
)
import scipy._lib.array_api_compat.numpy as np_compat

skip_xp_backends = pytest.mark.skip_xp_backends


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
        x, y = _asarray([0, 1, 2], xp=xp), _asarray(np.arange(3), xp=xp)
        ref = xp.asarray([0, 1, 2])
        xp_assert_equal(x, ref)
        xp_assert_equal(y, ref)

    @pytest.mark.filterwarnings("ignore: the matrix subclass")
    def test_raises(self):
        msg = "of type `numpy.ma.MaskedArray` are not supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace(np.ma.array(1), np.array(1))

        msg = "of type `numpy.matrix` are not supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace(np.array(1), np.matrix(1))

        msg = "only boolean and numerical dtypes are supported"
        with pytest.raises(TypeError, match=msg):
            array_namespace([object()])
        with pytest.raises(TypeError, match=msg):
            array_namespace('abc')

    def test_array_likes(self):
        # should be no exceptions
        array_namespace([0, 1, 2])
        array_namespace(1, 2, 3)
        array_namespace(1)

    @skip_xp_backends('jax.numpy',
                      reasons=["JAX arrays do not support item assignment"])
    @pytest.mark.usefixtures("skip_xp_backends")
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

    @array_api_compatible
    @pytest.mark.parametrize('dtype', ['int32', 'int64', 'float32', 'float64'])
    @pytest.mark.parametrize('shape', [(), (3,)])
    def test_strict_checks(self, xp, dtype, shape):
        # Check that `_strict_check` behaves as expected
        dtype = getattr(xp, dtype)
        x = xp.broadcast_to(xp.asarray(1, dtype=dtype), shape)
        x = x if shape else x[()]
        y = np_compat.asarray(1)[()]

        options = dict(check_namespace=True, check_dtype=False, check_shape=False)
        if xp == np:
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="Namespaces do not match."):
                xp_assert_equal(x, y, **options)

        options = dict(check_namespace=False, check_dtype=True, check_shape=False)
        if y.dtype.name in str(x.dtype):
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="dtypes do not match."):
                xp_assert_equal(x, y, **options)

        options = dict(check_namespace=False, check_dtype=False, check_shape=True)
        if x.shape == y.shape:
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="Shapes do not match."):
                xp_assert_equal(x, y, **options)

    @array_api_compatible
    def test_check_scalar(self, xp):
        if not is_numpy(xp):
            pytest.skip("Scalars only exist in NumPy")

        if is_numpy(xp):
            with pytest.raises(AssertionError, match="Types do not match."):
                xp_assert_equal(xp.asarray(0.), xp.float64(0))
            xp_assert_equal(xp.float64(0), xp.asarray(0.))


@array_api_compatible
class TestXP_Mean:
    @pytest.mark.parametrize('axis', [None, 1, -1, (-2, 2)])
    @pytest.mark.parametrize('weights', [None, True])
    @pytest.mark.parametrize('keepdims', [False, True])
    def test_xp_mean_basic(self, xp, axis, weights, keepdims):
        rng = np.random.default_rng(90359458245906)
        x = rng.random((3, 4, 5))
        x_xp = xp.asarray(x)
        w = w_xp = None

        if weights:
            w = rng.random((1, 5))
            w_xp = xp.asarray(w)
            x, w = np.broadcast_arrays(x, w)

        res = xp_mean(x_xp, weights=w_xp, axis=axis, keepdims=keepdims)
        ref = np.average(x, weights=w, axis=axis, keepdims=keepdims)

        xp_assert_close(res, xp.asarray(ref))

    def test_special_cases(self, xp):
        # non-broadcastable x and weights
        message = "...mismatch: objects cannot be broadcast to a single shape"
        x, w = xp.arange(10.), xp.zeros(5)
        with pytest.raises((ValueError, RuntimeError), match=message):
            xp_mean(x, weights=w)

        # weights sum to zero
        weights = xp.asarray([-1., 0., 1.])

        res = xp_mean(xp.asarray([1., 1., 1.]), weights=weights)
        xp_assert_close(res, xp.asarray(xp.nan))

        res = xp_mean(xp.asarray([2., 1., 1.]), weights=weights)
        xp_assert_close(res, xp.asarray(-np.inf))

        res = xp_mean(xp.asarray([1., 1., 2.]), weights=weights)
        xp_assert_close(res, xp.asarray(np.inf))

    def test_nan_policy(self, xp):
        x = xp.arange(10.)
        mask = (x == 3)
        x = xp.where(mask, xp.asarray(xp.nan), x)

        # nan_policy='raise' raises an error
        message = 'The input contains nan values'
        with pytest.raises(ValueError, match=message):
            xp_mean(x, nan_policy='raise')

        # `nan_policy='propagate'` is the default, and the result is NaN
        res1 = xp_mean(x)
        res2 = xp_mean(x, nan_policy='propagate')
        ref = xp.asarray(xp.nan)
        xp_assert_equal(res1, ref)
        xp_assert_equal(res2, ref)

        # `nan_policy='omit'` omits NaNs in `x`
        res = xp_mean(x, nan_policy='omit')
        ref = xp.mean(x[~mask])
        xp_assert_close(res, ref)

        # `nan_policy='omit'` omits NaNs in `weights`, too
        weights = xp.ones(10)
        weights = xp.where(mask, xp.asarray(xp.nan), weights)
        res = xp_mean(xp.arange(10.), weights=weights, nan_policy='omit')
        ref = xp.mean(x[~mask])
        xp_assert_close(res, ref)

        # Check for warning if omitting NaNs causes empty slice
        message = 'After omitting NaNs...'
        with pytest.warns(UserWarning, match=message):
            res = xp_mean(x * np.nan,  nan_policy='omit')
            ref = xp.asarray(xp.nan)
            xp_assert_equal(res, ref)

    def test_empty(self, xp):
        message = 'At least one slice along `axis` has...'

        with pytest.warns(UserWarning, match=message):
            res = xp_mean(xp.asarray([]))
            ref = xp.asarray(xp.nan)
            xp_assert_equal(res, ref)

        with pytest.warns(UserWarning, match=message):
            res = xp_mean(xp.asarray([[]]), axis=1)
            ref = xp.asarray([xp.nan])
            xp_assert_equal(res, ref)

        res = xp_mean(xp.asarray([[]]), axis=0)
        ref = xp.asarray([])
        xp_assert_equal(res, ref)

    def test_dtype(self, xp):
        max = xp.finfo(xp.float32).max
        x_np = np.asarray([max, max], dtype=np.float32)
        x_xp = xp.asarray(x_np)

        # Overflow occurs for float32 input
        with np.errstate(over='ignore'):
            res = xp_mean(x_xp)
            ref = np.mean(x_np)
            np.testing.assert_equal(ref, np.inf)
            xp_assert_close(res, xp.asarray(ref))

        # correct result is returned if `float64` is used
        res = xp_mean(x_xp, dtype=xp.float64)
        ref = xp.asarray(np.mean(np.asarray(x_np, dtype=np.float64)))
        xp_assert_close(res, ref)

    def test_integer(self, xp):
        # integer inputs are converted to the appropriate float
        x = xp.arange(10)
        y = xp.arange(10.)
        xp_assert_equal(xp_mean(x), xp_mean(y))
        xp_assert_equal(xp_mean(y, weights=x), xp_mean(y, weights=y))
