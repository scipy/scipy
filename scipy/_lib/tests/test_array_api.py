import re

import numpy as np
import pytest

from scipy._lib._array_api import (
    SCIPY_ARRAY_API, array_namespace, _asarray, xp_copy, xp_assert_equal, is_numpy,
    np_compat, xp_device, xp_promote, xp_result_device,
    xp_result_type, is_torch, _xp_copy_to_numpy
)
from scipy._external import array_api_extra as xpx
from scipy._lib._array_api_no_0d import xp_assert_equal as xp_assert_equal_no_0d
from scipy._external.array_api_extra.testing import lazy_xp_function

# Run all tests in this module in the Array API CI,
# including those without the `xp` fixture
pytestmark = pytest.mark.array_api_backends

lazy_xp_function(_asarray)
lazy_xp_function(xp_copy)


@pytest.mark.uses_xp_capabilities(False, reason="not applicable")
@pytest.mark.skipif(not SCIPY_ARRAY_API,
        reason="Array API test; set environment variable SCIPY_ARRAY_API=1 to run it")
class TestArrayAPI:

    def test_array_namespace(self):
        x, y = np.array([0, 1, 2]), np.array([0, 1, 2])
        xp = array_namespace(x, y)
        assert 'array_api_compat.numpy' in xp.__name__

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

    @pytest.mark.skip_xp_backends(np_only=True, reason="Array-likes")
    def test_array_likes(self, xp):
        """Test that if all parameters of array_namespace are Array-likes,
        the output is array_api_compat.numpy
        """
        assert array_namespace([0, 1, 2]) is xp
        assert array_namespace((0, 1, 2)) is xp
        assert array_namespace(1, 2, 3) is xp
        assert array_namespace(1) is xp
        assert array_namespace(np.int64(1)) is xp
        assert array_namespace([0, 1, 2], 3) is xp
        assert array_namespace() is xp
        assert array_namespace(None) is xp
        assert array_namespace(1, None) is xp
        assert array_namespace(None, 1) is xp

        # This only works when xp is numpy!
        assert array_namespace(np.asarray([1, 2]), [3, 4]) is xp
        assert array_namespace(np.int64(1), [3, 4]) is xp

    def test_array_and_array_likes_mix(self, xp):
        """Test that if there is at least one Array API object among
        the parameters of array_namespace, and all other parameters
        are scalars, the output is its namespace.

        If there are non-scalar Array-Likes, raise as in array-api-compat.
        """
        x = xp.asarray(1)
        assert array_namespace(x) is xp
        assert array_namespace(x, 1) is xp
        assert array_namespace(1, x) is xp
        assert array_namespace(None, x) is xp

        if is_numpy(xp):
            assert array_namespace(x, [1, 2]) is xp
        else:
            with pytest.raises(TypeError, match="Multiple namespaces"):
                array_namespace(x, [1, 2])
            with pytest.raises(TypeError, match="Multiple namespaces"):
                array_namespace(x, np.int64(1))
            with pytest.raises(TypeError, match="Multiple namespaces"):
                 # Subclass of float; matches array_api_compat behavior
                array_namespace(x, np.float64(1))
            with pytest.raises(TypeError, match="Multiple namespaces"):
                # Subclass of complex; matches array_api_compat behavior
                array_namespace(x, np.complex128(1))

    def test_array_api_extra_hook(self):
        """Test that the `array_namespace` function used by
        array-api-extra has been overridden by scipy
        """
        msg = "only boolean and numerical dtypes are supported"
        with pytest.raises(TypeError, match=msg):
            xpx.atleast_nd("abc", ndim=0)

    def test_jax_zero_gradient_array(self):
        """Test array_namespace special case for JAX zero-gradient arrays, which are
        numpy arrays but must be treated as JAX arrays.
        See matching code and tests in array_api_compat.
        """
        jax = pytest.importorskip("jax")
        xp = pytest.importorskip("jax.numpy")
        # Create numpy array with dtype=jax.float0
        jax_zero = jax.vmap(jax.grad(xp.float32, allow_int=True))(xp.arange(4))
        assert array_namespace(jax_zero) is xp

    def test_void_but_not_jax_zero_gradient_array(self):
        """A void dtype that is not a jax.float0 must not be caught in the
        special case for JAX zero-gradient arrays.
        """
        void = np.empty(0, dtype=np.dtype([]))
        with pytest.raises(TypeError, match="only boolean and numerical dtypes"):
            array_namespace(void)
        with pytest.raises(TypeError, match="only boolean and numerical dtypes"):
            array_namespace([void, void])

    def test_copy(self, xp):
        for _xp in [xp, None]:
            x = xp.asarray([1, 2, 3])
            y = xp_copy(x, xp=_xp)
            # with numpy we'd want to use np.shared_memory, but that's not specified
            # in the array-api
            assert id(x) != id(y)
            try:
                y[0] = 10
            except (TypeError, ValueError):
                pass
            else:
                assert x[0] != y[0]

    @pytest.mark.parametrize(
        "dtype",
        ["float32", "float64", "complex64", "complex128", "int32", "int64"],
    )
    @pytest.mark.parametrize(
        "data", [[], 1, [1, 2, 3], [[1, 2], [2, 3]]],
    )
    def test_copy_to_numpy(self, xp, data, dtype):
        xp_dtype = getattr(xp, dtype)
        np_dtype = getattr(np, dtype)
        x = xp.asarray(data, dtype=xp_dtype)
        y = _xp_copy_to_numpy(x)
        assert isinstance(y, np.ndarray)
        assert y.dtype == np_dtype
        assert x.shape == y.shape
        np.testing.assert_equal(y, np.asarray(data, dtype=np_dtype))
        if is_numpy(xp):
            # Ensure y is a copy when xp is numpy.
            assert id(x) != id(y)


    @pytest.mark.parametrize('dtype', ['int32', 'int64', 'float32', 'float64'])
    @pytest.mark.parametrize('shape', [(), (3,)])
    def test_strict_checks(self, xp, dtype, shape):
        # Check that `_strict_check` behaves as expected
        dtype = getattr(xp, dtype)
        x = xp.broadcast_to(xp.asarray(1, dtype=dtype), shape)
        x = x if shape else x[()]
        y = np_compat.asarray(1)[()]

        kwarg_names = ["check_namespace", "check_dtype", "check_shape", "check_0d"]
        options = dict(zip(kwarg_names, [True, False, False, False]))
        if is_numpy(xp):
            xp_assert_equal(x, y, **options)
        else:
            # `desired` can be of a different namespace
            # as long as `actual` matches the expectation set by `default_xp`
            xp_assert_equal(x, y, **options)
            with pytest.raises(
                AssertionError,
                match="Input does not have the desired array namespace",
            ):
                xp_assert_equal(y, x, **options)

        options = dict(zip(kwarg_names, [False, True, False, False]))
        if y.dtype.name in str(x.dtype):
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="dtypes do not match."):
                xp_assert_equal(x, y, **options)

        options = dict(zip(kwarg_names, [False, False, True, False]))
        if x.shape == y.shape:
            xp_assert_equal(x, y, **options)
        else:
            with pytest.raises(AssertionError, match="Shapes do not match."):
                xp_assert_equal(x, xp.asarray(y), **options)

        options = dict(zip(kwarg_names, [False, False, False, True]))
        if is_numpy(xp) and x.shape == y.shape:
            xp_assert_equal(x, y, **options)
        elif is_numpy(xp):
            with pytest.raises(AssertionError, match="Array-ness does not match."):
                xp_assert_equal(x, y, **options)

    @pytest.mark.skip_xp_backends(np_only=True, reason="Scalars only exist in NumPy")
    def test_check_scalar(self, xp):
        # identity always passes
        xp_assert_equal(xp.float64(0), xp.float64(0))
        xp_assert_equal(xp.asarray(0.), xp.asarray(0.))
        xp_assert_equal(xp.float64(0), xp.float64(0), check_0d=False)
        xp_assert_equal(xp.asarray(0.), xp.asarray(0.), check_0d=False)

        # Check default convention: 0d-arrays are distinguished from scalars
        message = "Array-ness does not match:.*"
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal(xp.asarray(0.), xp.float64(0))
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal(xp.float64(0), xp.asarray(0.))
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal(xp.asarray(42), xp.int64(42))
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal(xp.int64(42), xp.asarray(42))

        # with `check_0d=False`, scalars-vs-0d passes (if values match)
        xp_assert_equal(xp.asarray(0.), xp.float64(0), check_0d=False)
        xp_assert_equal(xp.float64(0), xp.asarray(0.), check_0d=False)
        # also with regular python objects
        xp_assert_equal(xp.asarray(0.), 0., check_0d=False)
        xp_assert_equal(0., xp.asarray(0.), check_0d=False)
        xp_assert_equal(xp.asarray(42), 42, check_0d=False)
        xp_assert_equal(42, xp.asarray(42), check_0d=False)

        # as an alternative to `check_0d=False`, explicitly expect scalar
        xp_assert_equal(xp.float64(0), xp.asarray(0.)[()])

    @pytest.mark.skip_xp_backends(np_only=True, reason="Scalars only exist in NumPy")
    def test_check_scalar_no_0d(self, xp):
        # identity passes, if first argument is not 0d (or check_0d=True)
        xp_assert_equal_no_0d(xp.float64(0), xp.float64(0))
        xp_assert_equal_no_0d(xp.float64(0), xp.float64(0), check_0d=True)
        xp_assert_equal_no_0d(xp.asarray(0.), xp.asarray(0.), check_0d=True)

        # by default, 0d values are forbidden as the first argument
        message = "Result is a NumPy 0d-array.*"
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal_no_0d(xp.asarray(0.), xp.asarray(0.))
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal_no_0d(xp.asarray(0.), xp.float64(0))
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal_no_0d(xp.asarray(42), xp.int64(42))

        # Check default convention: 0d-arrays are NOT distinguished from scalars
        xp_assert_equal_no_0d(xp.float64(0), xp.asarray(0.))
        xp_assert_equal_no_0d(xp.int64(42), xp.asarray(42))

        # opt in to 0d-check remains possible
        message = "Array-ness does not match:.*"
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal_no_0d(xp.asarray(0.), xp.float64(0), check_0d=True)
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal_no_0d(xp.float64(0), xp.asarray(0.), check_0d=True)
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal_no_0d(xp.asarray(42), xp.int64(0), check_0d=True)
        with pytest.raises(AssertionError, match=message):
            xp_assert_equal_no_0d(xp.int64(0), xp.asarray(42), check_0d=True)

        # scalars-vs-0d passes (if values match) also with regular python objects
        xp_assert_equal_no_0d(0., xp.asarray(0.))
        xp_assert_equal_no_0d(42, xp.asarray(42))


scalars = [1, 1., 1. + 1j]
lists = [[1], [1.], [1. + 1j]]
types = ('int8 int16 int32 int64 '
         'uint8 uint16 uint32 uint64 '
         'float32 float64 complex64 complex128').split()
arrays = [np.asarray([1], dtype=getattr(np, t)) for t in types]


def convert_type(x, xp):
    # Convert NumPy array to xp-array
    # Convert string to indicated dtype from xp
    # Return Python scalars unchanged
    if isinstance(x, np.ndarray):
        return xp.asarray(x)
    elif isinstance(x, str):
        return getattr(xp, x)
    return x


def is_inexact(x, xp):
    # Determine whether `x` is of inexact (real of complex floating) dtype
    x = xp.asarray(x) if np.isscalar(x) or isinstance(x, list) else x
    dtype = getattr(x, 'dtype', x)
    return xp.isdtype(dtype, ('real floating', 'complex floating'))


@pytest.mark.parametrize('x', scalars + lists + types + arrays)
@pytest.mark.parametrize('y', scalars + lists + types + arrays)
@pytest.mark.uses_xp_capabilities(False, reason="not applicable")
def test_xp_result_type_no_force(x, y, xp):
    # When force_floating==False (default), behavior of `xp_result_type`
    # should match that of `xp.result_type` on the same arguments after
    # converting lists to arrays of type `xp`.
    x = convert_type(x, xp)
    y = convert_type(y, xp)
    x_ref = xp.asarray(x) if isinstance(x, list) else x
    y_ref = xp.asarray(y) if isinstance(y, list) else y

    try:
        dtype_ref = xp.result_type(x_ref, y_ref)
        expected_error = None
    except Exception as e:
        expected_error = (type(e), str(e))

    if expected_error is not None:
        with pytest.raises(expected_error[0], match=re.escape(expected_error[1])):
            xp_result_type(x, y, xp=xp)
        return

    dtype_res = xp_result_type(x, y, xp=xp)
    assert dtype_res == dtype_ref


@pytest.mark.parametrize('x', scalars + lists + types + arrays)
@pytest.mark.parametrize('y', scalars + lists + types + arrays)
@pytest.mark.uses_xp_capabilities(False, reason="not applicable")
def test_xp_result_type_force_floating(x, y, xp):
    # When `force_floating==True`, behavior of `xp_result_type`
    # should match that of `xp.result_type` with `1.0` appended to the set of
    # arguments (after converting lists to arrays of type `xp`).
    # If this raises a `TypeError`, which is the case when the result
    # type is not defined by the standard, the result type should be
    # the result type of any inexact (real or complex floating) arguments
    # and the default floating point type.
    if (is_torch(xp) and not(isinstance(x, str) or isinstance(y, str))
            and np.isscalar(x) and np.isscalar(y)):
        pytest.skip("See 3/27/2024 comment at  data-apis/array-api-compat#277")

    x = convert_type(x, xp)
    y = convert_type(y, xp)
    x_ref = xp.asarray(x) if isinstance(x, list) else x
    y_ref = xp.asarray(y) if isinstance(y, list) else y

    expected_error = None
    try:
        dtype_ref = xp.result_type(x_ref, y_ref, 1.0)
    except TypeError:
        args = []
        if is_inexact(x_ref, xp):
            args.append(x_ref)
        if is_inexact(y_ref, xp):
            args.append(y_ref)
        dtype_ref = xp.result_type(*args, xp.asarray(1.0))
    except Exception as e:
        expected_error = (type(e), str(e))

    if expected_error is not None:
        with pytest.raises(expected_error[0], match=expected_error[1]):
            xp_result_type(x, y, xp=xp)
        return

    dtype_res = xp_result_type(x, y, force_floating=True, xp=xp)
    assert dtype_res == dtype_ref


@pytest.mark.uses_xp_capabilities(False, reason="not applicable")
def test_xp_promote_device(xp, devices):
    # Scalars and array-like iterables promoted alongside an array are
    # created on that array's device; array arguments keep their own device.
    # See gh-22680.
    for device in devices:
        x = xp.asarray([1., 2., 3.], device=device)

        # scalar rides along on the array's device
        x2, scalar = xp_promote(x, 1.5, xp=xp)
        assert xp_device(x2) == xp_device(x)
        assert xp_device(scalar) == xp_device(x)

        # so does a list, and broadcasting preserves the device
        x2, lst = xp_promote(x, [1., 2., 3.], broadcast=True, xp=xp)
        assert xp_device(lst) == xp_device(x)

        # order does not matter: device comes from the (first) array argument
        scalar, x2 = xp_promote(2, x, force_floating=True, xp=xp)
        assert xp_device(scalar) == xp_device(x)

        # with arrays on several devices, the *first* array argument is the
        # anchor for scalars; each array keeps its own device
        y = xp.asarray([4., 5., 6.])  # on the default device
        scalar, x2, y2 = xp_promote(2, x, y, force_floating=True, xp=xp)
        assert xp_device(scalar) == xp_device(x)
        assert xp_device(x2) == xp_device(x)
        assert xp_device(y2) == xp_device(y)

        # None passes through; remaining args still promoted on device
        x2, none, scalar = xp_promote(x, None, 0.5, xp=xp)
        assert none is None
        assert xp_device(scalar) == xp_device(x)

        # an array argument keeps its own device
        x2 = xp_promote(x, force_floating=True, xp=xp)
        assert xp_device(x2) == xp_device(x)

        # NumPy scalars and arrays report device 'cpu' but are *host data*:
        # they must ride along on the array argument's device like python
        # scalars do (they land on the default device otherwise). This mix
        # is user-reachable: `array_namespace` rejects NumPy *sample* arrays
        # mixed with another backend's arrays, but NumPy scalars arrive as
        # numeric keyword parameters that are not validated as arrays, e.g.
        # `differentiate(f, x_gpu, initial_step=1/np.sqrt(2))`.
        x2, np_scalar = xp_promote(x, np.float64(0.5), xp=xp)
        assert xp_device(np_scalar) == xp_device(x)
        x2, np_arr = xp_promote(x, np.asarray([1.0, 2.0, 3.0]), xp=xp)
        assert xp_device(np_arr) == xp_device(x)


@pytest.mark.uses_xp_capabilities(False, reason="not applicable")
def test_xp_result_device(xp, devices):
    # `xp_result_device` anchors the result device on the first argument
    # that carries a device, skipping python scalars and NumPy host data.
    # No device-carrying argument -> None (create on the backend default).
    # NumPy arrays, and arrays of backends without a `.device` attribute
    # (e.g. Dask), are treated like host data: always None.
    for device in devices:
        x = xp.asarray([1., 2.], device=device)
        host_like = is_numpy(xp) or not hasattr(x, "device")
        expected = None if host_like else xp_device(x)
        assert xp_result_device(x) == expected
        assert xp_result_device(1.5, None, x) == expected
        assert xp_result_device(np.float64(0.5), x) == expected
        assert xp_result_device(np.asarray([1.0]), x) == expected

        # the *first* device-carrying argument is the anchor
        y = xp.asarray([3., 4.])  # on the default device
        expected_y = None if host_like else xp_device(y)
        assert xp_result_device(x, y) == expected
        assert xp_result_device(y, x) == expected_y
        assert xp_result_device(1.5, y, x) == expected_y
    assert xp_result_device() is None
    assert xp_result_device(1.5, None) is None
    assert xp_result_device(np.asarray([1.0])) is None
