"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""
import os

import numpy as np
from numpy.core.numerictypes import typecodes
# probably want to vendor it (submodule)
import array_api_compat
import array_api_compat.numpy

__all__ = ['array_namespace', 'as_xparray', 'as_xparray_namespace']


# SCIPY_ARRAY_API, array_api_dispatch is used by sklearn
array_api_dispatch = os.environ.get("array_api_dispatch", False)
SCIPY_ARRAY_API = os.environ.get("SCIPY_ARRAY_API", array_api_dispatch)

_GLOBAL_CONFIG = {"SCIPY_ARRAY_API": SCIPY_ARRAY_API}


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
            raise ValueError('object arrays are not supported')


def _check_finite(array, xp):
    """Check for NaNs or Infs."""
    # same as np.asarray_chkfinite
    if array.dtype.char in typecodes['AllFloat'] and not xp.isfinite(array).all():
        raise ValueError(
            "array must not contain infs or NaNs"
        )


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
        return np

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


def as_xparray_namespace(*arrays):
    """Validate and convert arrays to a common namespace.

    Parameters
    ----------
    *arrays : sequence of array_like
        Arrays to validate and convert.

    Returns
    -------
    *arrays : sequence of array_like
        Validated and converted arrays to the common namespace.
    namespace : module
        Common namespace.

    Notes
    -----
    This function is meant to be called from each public function in a SciPy
    submodule it does the following:

    1. Check for the global switch: SCIPY_ARRAY_API. This can also be accessed
       dynamically through ``_GLOBAL_CONFIG['SCIPY_ARRAY_API']``.
    2. `compliance_scipy` raise exceptions on known-bad subclasses. See
       it's definition for more details.
    3. Determine the namespace, without doing any coercion of array(-like)
       inputs.
    4. Call `xp.asarray` on all array.

    Examples
    --------
    >>> import numpy as np
    >>> x, y, xp = as_xparray_namespace(np.array([0, 1, 2]), np.array([0, 1, 2]))
    >>> xp.__name__
    'array_api_compat.numpy'
    >>> x, y
    (array([0, 1, 2]), array([0, 1, 2]))

    """
    arrays = list(arrays)
    xp = array_namespace(*arrays)

    for i, array in enumerate(arrays):
        arrays[i] = xp.asarray(array)

    return *arrays, xp


def to_numpy(array, xp):
    """Convert `array` into a NumPy ndarray on the CPU.

    ONLY FOR TESTING
    """
    xp_name = xp.__name__

    if xp_name in {"array_api_compat.torch", "torch"}:
        return array.cpu().numpy()
    elif xp_name == "cupy.array_api":
        return array._array.get()
    elif xp_name in {"array_api_compat.cupy", "cupy"}:  # pragma: nocover
        return array.get()

    return np.asarray(array)
