"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""
from __future__ import annotations

import os

import numpy as np
import scipy._lib.array_api_compat.array_api_compat as array_api_compat
from scipy._lib.array_api_compat.array_api_compat import size
import scipy._lib.array_api_compat.array_api_compat.numpy as array_api_compat_numpy

__all__ = ['array_namespace', 'as_xparray', 'size']


# SCIPY_ARRAY_API, array_api_dispatch is used by sklearn
array_api_dispatch = os.environ.get("array_api_dispatch", False)
SCIPY_ARRAY_API: str | bool = os.environ.get("SCIPY_ARRAY_API", array_api_dispatch)
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
    array, dtype=None, order=None, copy=None, *, xp=None, check_finite=True
):
    """Drop-in replacement for `np.asarray`.

    Memory layout parameter `order` is not exposed in the Array API standard.
    `order` is only enforced if the input array implementation
    is NumPy based, otherwise `order` is just silently ignored.

    The purpose of this helper is to make it possible to share code for data
    container validation without memory copies for both downstream use cases.
    """
    if xp is None:
        xp = array_namespace(array)
    if xp.__name__ in {"numpy", "array_api_compat.numpy", "numpy.array_api"}:
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
