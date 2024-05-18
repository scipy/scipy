"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""
from __future__ import annotations

import os
import warnings

from types import ModuleType
from typing import Any, Literal, TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from scipy._lib import array_api_compat
from scipy._lib.array_api_compat import (
    is_array_api_obj,
    size,
    numpy as np_compat,
)

__all__ = ['array_namespace', '_asarray', 'size']


# To enable array API and strict array-like input validation
SCIPY_ARRAY_API: str | bool = os.environ.get("SCIPY_ARRAY_API", False)
# To control the default device - for use in the test suite only
SCIPY_DEVICE = os.environ.get("SCIPY_DEVICE", "cpu")

_GLOBAL_CONFIG = {
    "SCIPY_ARRAY_API": SCIPY_ARRAY_API,
    "SCIPY_DEVICE": SCIPY_DEVICE,
}


if TYPE_CHECKING:
    Array = Any  # To be changed to a Protocol later (see array-api#589)
    ArrayLike = Array | npt.ArrayLike


def compliance_scipy(arrays: list[ArrayLike]) -> list[Array]:
    """Raise exceptions on known-bad subclasses.

    The following subclasses are not supported and raise and error:
    - `numpy.ma.MaskedArray`
    - `numpy.matrix`
    - NumPy arrays which do not have a boolean or numerical dtype
    - Any array-like which is neither array API compatible nor coercible by NumPy
    - Any array-like which is coerced by NumPy to an unsupported dtype
    """
    for i in range(len(arrays)):
        array = arrays[i]
        if isinstance(array, np.ma.MaskedArray):
            raise TypeError("Inputs of type `numpy.ma.MaskedArray` are not supported.")
        elif isinstance(array, np.matrix):
            raise TypeError("Inputs of type `numpy.matrix` are not supported.")
        if isinstance(array, (np.ndarray, np.generic)):
            dtype = array.dtype
            if not (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bool_)):
                raise TypeError(f"An argument has dtype `{dtype!r}`; "
                                f"only boolean and numerical dtypes are supported.")
        elif not is_array_api_obj(array):
            try:
                array = np.asanyarray(array)
            except TypeError:
                raise TypeError("An argument is neither array API compatible nor "
                                "coercible by NumPy.")
            dtype = array.dtype
            if not (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bool_)):
                message = (
                    f"An argument was coerced to an unsupported dtype `{dtype!r}`; "
                    f"only boolean and numerical dtypes are supported."
                )
                raise TypeError(message)
            arrays[i] = array
    return arrays


def _check_finite(array: Array, xp: ModuleType) -> None:
    """Check for NaNs or Infs."""
    msg = "array must not contain infs or NaNs"
    try:
        if not xp.all(xp.isfinite(array)):
            raise ValueError(msg)
    except TypeError:
        raise ValueError(msg)


def array_namespace(*arrays: Array) -> ModuleType:
    """Get the array API compatible namespace for the arrays xs.

    Parameters
    ----------
    *arrays : sequence of array_like
        Arrays used to infer the common namespace.

    Returns
    -------
    namespace : module
        Common namespace.

    Notes
    -----
    Thin wrapper around `array_api_compat.array_namespace`.

    1. Check for the global switch: SCIPY_ARRAY_API. This can also be accessed
       dynamically through ``_GLOBAL_CONFIG['SCIPY_ARRAY_API']``.
    2. `compliance_scipy` raise exceptions on known-bad subclasses. See
       its definition for more details.

    When the global switch is False, it defaults to the `numpy` namespace.
    In that case, there is no compliance check. This is a convenience to
    ease the adoption. Otherwise, arrays must comply with the new rules.
    """
    if not _GLOBAL_CONFIG["SCIPY_ARRAY_API"]:
        # here we could wrap the namespace if needed
        return np_compat

    _arrays = [array for array in arrays if array is not None]

    _arrays = compliance_scipy(_arrays)

    return array_api_compat.array_namespace(*_arrays)


def _asarray(
        array: ArrayLike,
        dtype: Any = None,
        order: Literal['K', 'A', 'C', 'F'] | None = None,
        copy: bool | None = None,
        *,
        xp: ModuleType | None = None,
        check_finite: bool = False,
        subok: bool = False,
    ) -> Array:
    """SciPy-specific replacement for `np.asarray` with `order`, `check_finite`, and
    `subok`.

    Memory layout parameter `order` is not exposed in the Array API standard.
    `order` is only enforced if the input array implementation
    is NumPy based, otherwise `order` is just silently ignored.

    `check_finite` is also not a keyword in the array API standard; included
    here for convenience rather than that having to be a separate function
    call inside SciPy functions.

    `subok` is included to allow this function to preserve the behaviour of
    `np.asanyarray` for NumPy based inputs.
    """
    if xp is None:
        xp = array_namespace(array)
    if xp.__name__ in {"numpy", "scipy._lib.array_api_compat.numpy"}:
        # Use NumPy API to support order
        if copy is True:
            array = np.array(array, order=order, dtype=dtype, subok=subok)
        elif subok:
            array = np.asanyarray(array, order=order, dtype=dtype)
        else:
            array = np.asarray(array, order=order, dtype=dtype)

        # At this point array is a NumPy ndarray. We convert it to an array
        # container that is consistent with the input's namespace.
        array = xp.asarray(array)
    else:
        try:
            array = xp.asarray(array, dtype=dtype, copy=copy)
        except TypeError:
            coerced_xp = array_namespace(xp.asarray(3))
            array = coerced_xp.asarray(array, dtype=dtype, copy=copy)

    if check_finite:
        _check_finite(array, xp)

    return array


def atleast_nd(x: Array, *, ndim: int, xp: ModuleType | None = None) -> Array:
    """Recursively expand the dimension to have at least `ndim`."""
    if xp is None:
        xp = array_namespace(x)
    x = xp.asarray(x)
    if x.ndim < ndim:
        x = xp.expand_dims(x, axis=0)
        x = atleast_nd(x, ndim=ndim, xp=xp)
    return x


def copy(x: Array, *, xp: ModuleType | None = None) -> Array:
    """
    Copies an array.

    Parameters
    ----------
    x : array

    xp : array_namespace

    Returns
    -------
    copy : array
        Copied array

    Notes
    -----
    This copy function does not offer all the semantics of `np.copy`, i.e. the
    `subok` and `order` keywords are not used.
    """
    # Note: xp.asarray fails if xp is numpy.
    if xp is None:
        xp = array_namespace(x)

    return _asarray(x, copy=True, xp=xp)


def is_numpy(xp: ModuleType) -> bool:
    return xp.__name__ in ('numpy', 'scipy._lib.array_api_compat.numpy')


def is_cupy(xp: ModuleType) -> bool:
    return xp.__name__ in ('cupy', 'scipy._lib.array_api_compat.cupy')


def is_torch(xp: ModuleType) -> bool:
    return xp.__name__ in ('torch', 'scipy._lib.array_api_compat.torch')

def is_jax(xp):
    return xp.__name__ in ('jax.numpy', 'jax.experimental.array_api')


def _strict_check(actual, desired, xp,
                  check_namespace=True, check_dtype=True, check_shape=True):
    __tracebackhide__ = True  # Hide traceback for py.test
    if check_namespace:
        _assert_matching_namespace(actual, desired)

    desired = xp.asarray(desired)

    if check_dtype:
        _msg = f"dtypes do not match.\nActual: {actual.dtype}\nDesired: {desired.dtype}"
        assert actual.dtype == desired.dtype, _msg

    if check_shape:
        _msg = f"Shapes do not match.\nActual: {actual.shape}\nDesired: {desired.shape}"
        assert actual.shape == desired.shape, _msg
        _check_scalar(actual, desired, xp)

    desired = xp.broadcast_to(desired, actual.shape)
    return desired


def _assert_matching_namespace(actual, desired):
    __tracebackhide__ = True  # Hide traceback for py.test
    actual = actual if isinstance(actual, tuple) else (actual,)
    desired_space = array_namespace(desired)
    for arr in actual:
        arr_space = array_namespace(arr)
        _msg = (f"Namespaces do not match.\n"
                f"Actual: {arr_space.__name__}\n"
                f"Desired: {desired_space.__name__}")
        assert arr_space == desired_space, _msg


def _check_scalar(actual, desired, xp):
    __tracebackhide__ = True  # Hide traceback for py.test
    # Shape check alone is sufficient unless desired.shape == (). Also,
    # only NumPy distinguishes between scalars and arrays.
    if desired.shape != () or not is_numpy(xp):
        return
    # We want to follow the conventions of the `xp` library. Libraries like
    # NumPy, for which `np.asarray(0)[()]` returns a scalar, tend to return
    # a scalar even when a 0D array might be more appropriate:
    # import numpy as np
    # np.mean([1, 2, 3])  # scalar, not 0d array
    # np.asarray(0)*2  # scalar, not 0d array
    # np.sin(np.asarray(0))  # scalar, not 0d array
    # Libraries like CuPy, for which `cp.asarray(0)[()]` returns a 0D array,
    # tend to return a 0D array in scenarios like those above.
    # Therefore, regardless of whether the developer provides a scalar or 0D
    # array for `desired`, we would typically want the type of `actual` to be
    # the type of `desired[()]`. If the developer wants to override this
    # behavior, they can set `check_shape=False`.
    desired = desired[()]
    _msg = f"Types do not match:\n Actual: {type(actual)}\n Desired: {type(desired)}"
    assert (xp.isscalar(actual) and xp.isscalar(desired)
            or (not xp.isscalar(actual) and not xp.isscalar(desired))), _msg


def xp_assert_equal(actual, desired, check_namespace=True, check_dtype=True,
                    check_shape=True, err_msg='', xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)
    desired = _strict_check(actual, desired, xp, check_namespace=check_namespace,
                            check_dtype=check_dtype, check_shape=check_shape)
    if is_cupy(xp):
        return xp.testing.assert_array_equal(actual, desired, err_msg=err_msg)
    elif is_torch(xp):
        # PyTorch recommends using `rtol=0, atol=0` like this
        # to test for exact equality
        err_msg = None if err_msg == '' else err_msg
        return xp.testing.assert_close(actual, desired, rtol=0, atol=0, equal_nan=True,
                                       check_dtype=False, msg=err_msg)
    # JAX uses `np.testing`
    return np.testing.assert_array_equal(actual, desired, err_msg=err_msg)


def xp_assert_close(actual, desired, rtol=None, atol=0, check_namespace=True,
                    check_dtype=True, check_shape=True, err_msg='', xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)
    desired = _strict_check(actual, desired, xp, check_namespace=check_namespace,
                            check_dtype=check_dtype, check_shape=check_shape)

    floating = xp.isdtype(actual.dtype, ('real floating', 'complex floating'))
    if rtol is None and floating:
        # multiplier of 4 is used as for `np.float64` this puts the default `rtol`
        # roughly half way between sqrt(eps) and the default for
        # `numpy.testing.assert_allclose`, 1e-7
        rtol = xp.finfo(actual.dtype).eps**0.5 * 4
    elif rtol is None:
        rtol = 1e-7

    if is_cupy(xp):
        return xp.testing.assert_allclose(actual, desired, rtol=rtol,
                                          atol=atol, err_msg=err_msg)
    elif is_torch(xp):
        err_msg = None if err_msg == '' else err_msg
        return xp.testing.assert_close(actual, desired, rtol=rtol, atol=atol,
                                       equal_nan=True, check_dtype=False, msg=err_msg)
    # JAX uses `np.testing`
    return np.testing.assert_allclose(actual, desired, rtol=rtol,
                                      atol=atol, err_msg=err_msg)


def xp_assert_less(actual, desired, check_namespace=True, check_dtype=True,
                   check_shape=True, err_msg='', verbose=True, xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)
    desired = _strict_check(actual, desired, xp, check_namespace=check_namespace,
                            check_dtype=check_dtype, check_shape=check_shape)
    if is_cupy(xp):
        return xp.testing.assert_array_less(actual, desired,
                                            err_msg=err_msg, verbose=verbose)
    elif is_torch(xp):
        if actual.device.type != 'cpu':
            actual = actual.cpu()
        if desired.device.type != 'cpu':
            desired = desired.cpu()
    # JAX uses `np.testing`
    return np.testing.assert_array_less(actual, desired,
                                        err_msg=err_msg, verbose=verbose)


def cov(x: Array, *, xp: ModuleType | None = None) -> Array:
    if xp is None:
        xp = array_namespace(x)

    X = copy(x, xp=xp)
    dtype = xp.result_type(X, xp.float64)

    X = atleast_nd(X, ndim=2, xp=xp)
    X = xp.asarray(X, dtype=dtype)

    avg = xp.mean(X, axis=1)
    fact = X.shape[1] - 1

    if fact <= 0:
        warnings.warn("Degrees of freedom <= 0 for slice",
                      RuntimeWarning, stacklevel=2)
        fact = 0.0

    X -= avg[:, None]
    X_T = X.T
    if xp.isdtype(X_T.dtype, 'complex floating'):
        X_T = xp.conj(X_T)
    c = X @ X_T
    c /= fact
    axes = tuple(axis for axis, length in enumerate(c.shape) if length == 1)
    return xp.squeeze(c, axis=axes)


def xp_unsupported_param_msg(param: Any) -> str:
    return f'Providing {param!r} is only supported for numpy arrays.'


def is_complex(x: Array, xp: ModuleType) -> bool:
    return xp.isdtype(x.dtype, 'complex floating')

def scipy_namespace_for(xp):
    """
    Return the `scipy` namespace for alternative backends, where it exists,
    such as `cupyx.scipy` and `jax.scipy`. Useful for ad hoc dispatching.

    Default: return `scipy` (this package).
    """


    if is_cupy(xp):
        import cupyx  # type: ignore[import-not-found]
        return cupyx.scipy

    if is_jax(xp):
        import jax  # type: ignore[import-not-found]
        return jax.scipy

    import scipy
    return scipy

# temporary substitute for xp.minimum, which is not yet in all backends
# or covered by array_api_compat.
def xp_minimum(x1, x2):
    # xp won't be passed in because it doesn't need to be passed in to xp.minimum
    xp = array_namespace(x1, x2)
    if hasattr(xp, 'minimum'):
        return xp.minimum(x1, x2)
    x1, x2 = xp.broadcast_arrays(x1, x2)
    dtype = xp.result_type(x1.dtype, x2.dtype)
    res = xp.asarray(x1, copy=True, dtype=dtype)
    i = (x2 < x1) | xp.isnan(x2)
    res[i] = x2[i]
    return res[()] if res.ndim == 0 else res


# temporary substitute for xp.clip, which is not yet in all backends
# or covered by array_api_compat.
def xp_clip(x, a, b, xp=None):
    xp = array_namespace(x) if xp is None else xp
    a, b = xp.asarray(a, dtype=x.dtype), xp.asarray(b, dtype=x.dtype)
    if hasattr(xp, 'clip'):
        return xp.clip(x, a, b)
    x, a, b = xp.broadcast_arrays(x, a, b)
    y = xp.asarray(x, copy=True)
    ia = y < a
    y[ia] = a[ia]
    ib = y > b
    y[ib] = b[ib]
    return y[()] if y.ndim == 0 else y


# temporary substitute for xp.moveaxis, which is not yet in all backends
# or covered by array_api_compat.
def xp_moveaxis_to_end(x, source, xp=None):
    xp = array_namespace(xp) if xp is None else xp
    axes = list(range(x.ndim))
    temp = axes.pop(source)
    axes = axes + [temp]
    return xp.permute_dims(x, axes)


# temporary substitute for xp.copysign, which is not yet in all backends
# or covered by array_api_compat.
def xp_copysign(x1, x2, xp=None):
    # no attempt to account for special cases
    xp = array_namespace(x1, x2) if xp is None else xp
    abs_x1 = xp.abs(x1)
    return xp.where(x2 >= 0, abs_x1, -abs_x1)


# partial substitute for xp.sign, which does not cover the NaN special case
# that I need. (https://github.com/data-apis/array-api-compat/issues/136)
def xp_sign(x, xp=None):
    xp = array_namespace(x) if xp is None else xp
    if is_numpy(xp):  # only NumPy implements the special cases correctly
        return xp.sign(x)
    sign = xp.full_like(x, xp.nan)
    one = xp.asarray(1, dtype=x.dtype)
    sign = xp.where(x > 0, one, sign)
    sign = xp.where(x < 0, -one, sign)
    sign = xp.where(x == 0, 0*one, sign)
    return sign


def xp_add_reduced_axes(res, axis, initial_shape, *, xp=None):
    xp = array_namespace(res) if xp is None else xp

    if axis is None:
        final_shape = (1,) * len(initial_shape)
    else:
        # axis can be a scalar or sequence
        axes = (axis,) if xp.asarray(axis).ndim == 0 else axis
        final_shape = list(initial_shape)
        for i in axes:
            final_shape[i] = 1

    return xp.reshape(res, final_shape)


# array-API compatible substitute for np.mean, np.nanmean, np.average
def xp_mean(x, *, axis=None, weights=None, keepdims=False, nan_policy='propagate',
            dtype=None, xp=None):
    r"""Compute the arithmetic mean along the specified axis.

    Parameters
    ----------
    x : real floating array
        Array containing real numbers whose mean is desired.
    axis : int or tuple of ints, default: None
        If an int or tuple of ints, the axis or axes of the input along which
        to compute the statistic. The statistic of each axis-slice (e.g. row)
        of the input will appear in a corresponding element of the output.
        If ``None``, the input will be raveled before computing the statistic.
    weights : real floating array, optional
        If specified, an array of weights associated with the values in `x`;
        otherwise ``1``. If `weights` and `x` do not have the same shape, the
        arrays will be broadcasted before performing the calculation. See
        Notes for details.
    keepdims : boolean, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with length one. With this option,
        the result will broadcast correctly against the input array.
    nan_policy : {'propagate', 'omit', 'raise'}, default: 'propagate'
        Defines how to handle input NaNs.

        - ``propagate``: if a NaN is present in the axis slice (e.g. row) along
          which the  statistic is computed, the corresponding entry of the output
          will be NaN.
        - ``omit``: NaNs will be omitted when performing the calculation.
          If insufficient data remains in the axis slice along which the
          statistic is computed, the corresponding entry of the output will be
          NaN.
        - ``raise``: if a NaN is present, a ``ValueError`` will be raised.

    dtype : dtype, optional
        Type to use in computing the mean. For integer inputs, the default is
        the default float type of the array library; for floating point inputs,
        the dtype is that of the input.

    Returns
    -------
    out : array
        The mean of each slice

    Notes
    -----
    Let :math:`x_i` represent element :math:`i` of data `x` and let :math:`w_i`
    represent the corresponding element of `weights` after broadcasting. Then the
    mean is given by:

    .. math::

        \frac{ \sum_{i=0}^{n-1} w_i x_i }
             { \sum_{i=0}^{n-1} i w_i }

    where `n` is the number of elements along a slice. Note that this simplifies
    to the familiar :math:`(\sum_i x_i) / n` when the weights are all ``1`` (default).

    The behavior of this function with respect to weights is somewhat different
    from that of `np.average`. For instance,
    `np.average` raises an error when `axis` is not specified and the shapes of `x`
    and a `weights` array are not the same; `xp_mean` simply broadcasts the two.
    Also, `np.average` raises an error when weights sum to zero along a slice;
    `xp_mean` computes the appropriate result.

    Note that according to the formula, including NaNs with zero weights is not
    the same as *omitting* NaNs with `nan_policy='omit'`; in the former case,
    the NaNs will continue to propagate through the calculation whereas in the
    latter case, the NaNs are excluded entirely.

    """
    xp = array_namespace(x) if xp is None else xp
    x = xp.asarray(x, dtype=dtype)
    weights = xp.asarray(weights, dtype=dtype) if weights is not None else weights

    if not xp.isdtype(x.dtype, 'real floating'):
        dtype = xp.asarray(1.).dtype
        x = xp.asarray(x, dtype=dtype)
    if weights is not None and not xp.isdtype(weights.dtype, 'real floating'):
        dtype = xp.asarray(1.).dtype
        weights = xp.asarray(weights, dtype=dtype)

    if weights is not None and x.shape != weights.shape:
        x, weights = xp.broadcast_arrays(x, weights)

    message = ('At least one slice along `axis` has zero length; '
               'corresponding slices of the output will be NaN.')
    if size(x) == 0:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = xp.mean(x, axis=axis, keepdims=keepdims)
        if size(res) != 0:
            warnings.warn(message, UserWarning, stacklevel=2)
        return res

    # avoid circular import
    from scipy._lib._util import _contains_nan
    contains_nan, _ = _contains_nan(x, nan_policy, xp_ok=True, xp=xp)
    if weights is not None:
        contains_nan_w, _ = _contains_nan(weights, nan_policy, xp_ok=True, xp=xp)
        contains_nan = contains_nan | contains_nan_w

    message = ('After omitting NaNs, at least one slice along `axis` has zero '
               'length; corresponding slices of the output will be NaN.')
    if contains_nan and nan_policy == 'omit':
        i = xp.isnan(x)
        i = (i | xp.isnan(weights)) if weights is not None else i
        if xp.any(xp.all(i, axis=axis)):
            warnings.warn(message, UserWarning, stacklevel=2)
        weights = xp.ones_like(x) if weights is None else weights
        x = xp.where(i, xp.asarray(0, dtype=x.dtype), x)
        weights = xp.where(i, xp.asarray(0, dtype=x.dtype), weights)

    if weights is None:
        return xp.mean(x, axis=axis, keepdims=keepdims)

    norm = xp.sum(weights, axis=axis)
    wsum = xp.sum(x * weights, axis=axis)
    with np.errstate(divide='ignore', invalid='ignore'):
        res = wsum/norm

    res = xp_add_reduced_axes(res, axis, x.shape, xp=xp) if keepdims else res
    return res[()] if res.ndim == 0 else res
