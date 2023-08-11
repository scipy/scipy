"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import scipy._lib.array_api_compat.array_api_compat as array_api_compat
from scipy._lib.array_api_compat.array_api_compat import size
import scipy._lib.array_api_compat.array_api_compat.numpy as array_api_compat_numpy

__all__ = ['array_namespace', 'as_xparray', 'size']


# To enable array API and strict array-like input validation
SCIPY_ARRAY_API: str | bool = os.environ.get("SCIPY_ARRAY_API", False)
# To control the default device - for use in the test suite only
SCIPY_DEVICE = os.environ.get("SCIPY_DEVICE", "cpu")

_GLOBAL_CONFIG = {
    "SCIPY_ARRAY_API": SCIPY_ARRAY_API,
    "SCIPY_DEVICE": SCIPY_DEVICE,
}


def compliance_scipy(*arrays):
    """Raise exceptions on known-bad subclasses.

    The following subclasses are not supported and raise and error:
    - `np.ma.MaskedArray`
    - `numpy.matrix`
    - Any array-like which is not Array API compatible
    """
    for array in arrays:
        if isinstance(array, np.ma.MaskedArray):
            raise TypeError("'numpy.ma.MaskedArray' are not supported")
        elif isinstance(array, np.matrix):
            raise TypeError("'numpy.matrix' are not supported")
        elif not array_api_compat.is_array_api_obj(array):
            raise TypeError("Only support Array API compatible arrays")
        elif array.dtype is np.dtype('O'):
            raise TypeError('object arrays are not supported')


def _check_finite(array, xp):
    """Check for NaNs or Infs."""
    msg = "array must not contain infs or NaNs"
    try:
        if not xp.all(xp.isfinite(array)):
            raise ValueError(msg)
    except TypeError:
        raise ValueError(msg)


def array_namespace(*arrays):
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
       it's definition for more details.

    When the global switch is False, it defaults to the `numpy` namespace.
    In that case, there is no compliance check. This is a convenience to
    ease the adoption. Otherwise, arrays must comply with the new rules.
    """
    if not _GLOBAL_CONFIG["SCIPY_ARRAY_API"]:
        # here we could wrap the namespace if needed
        return array_api_compat_numpy

    arrays = [array for array in arrays if array is not None]

    compliance_scipy(*arrays)

    return array_api_compat.array_namespace(*arrays)


def as_xparray(
    array, dtype=None, order=None, copy=None, *, xp=None, check_finite=False
):
    """SciPy-specific replacement for `np.asarray` with `order` and `check_finite`.

    Memory layout parameter `order` is not exposed in the Array API standard.
    `order` is only enforced if the input array implementation
    is NumPy based, otherwise `order` is just silently ignored.

    `check_finite` is also not a keyword in the array API standard; included
    here for convenience rather than that having to be a separate function
    call inside SciPy functions.
    """
    if xp is None:
        xp = array_namespace(array)
    if xp.__name__ in {"numpy", "scipy._lib.array_api_compat.array_api_compat.numpy"}:
        # Use NumPy API to support order
        if copy is True:
            array = np.array(array, order=order, dtype=dtype)
        else:
            array = np.asarray(array, order=order, dtype=dtype)

        # At this point array is a NumPy ndarray. We convert it to an array
        # container that is consistent with the input's namespace.
        array = xp.asarray(array)
    else:
        array = xp.asarray(array, dtype=dtype, copy=copy)

    if check_finite:
        _check_finite(array, xp)

    return array


def atleast_nd(x, *, ndim, xp):
    """Recursively expand the dimension to have at least `ndim`."""
    x = xp.asarray(x)
    if x.ndim < ndim:
        x = xp.expand_dims(x, axis=0)
        x = atleast_nd(x, ndim=ndim, xp=xp)
    return x


def copy(x, *, xp=None):
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

    return as_xparray(x, copy=True, xp=xp)


def cov(x, *, xp=None):
    if xp is None:
        xp = array_namespace(x)

    X = copy(x, xp=xp)
    dtype = xp.result_type(X, xp.float64)

    X = atleast_nd(X, ndim=2, xp=xp)
    X = xp.asarray(X, dtype=dtype)

    avg = xp.mean(X, axis=1)
    w_sum = xp.asarray(size(X) / size(avg), dtype=avg.dtype)
    if w_sum.shape != avg.shape:
        w_sum = copy(xp.broadcast_to(w_sum, avg.shape), xp=xp)

    w_sum = w_sum[0]

    fact = X.shape[1] - 1

    if fact <= 0:
        warnings.warn("Degrees of freedom <= 0 for slice",
                      RuntimeWarning, stacklevel=2)
        fact = 0.0

    X -= avg[:, None]
    X_T = X.T
    if xp.isdtype(X_T.dtype, 'complex floating'):
        X_T = xp.conj(X_T)
    c = xp.matmul(X, X_T)
    c *= 1 / fact
    axis = tuple(index for index, value in enumerate(c.shape) if value == 1)
    return xp.squeeze(c, axis=axis)
