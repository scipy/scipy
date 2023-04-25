"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""
import os

import numpy as np
# probably want to vendor it (submodule)
import array_api_compat
import array_api_compat.numpy

# SCIPY_ARRAY_API, array_api_dispatch is used by sklearn
USE_ARRAY_API = os.environ.get("array_api_dispatch", False)

__all__ = ['namespace_from_arrays', 'asarray', 'asarray_namespace']


def namespace_from_arrays(*arrays, single_namespace=True):
    # if we cannot get the namespace, np is used
    namespaces = set()
    for array in arrays:
        try:
            namespaces.add(array_api_compat.array_namespace(array))
        except TypeError:
            namespaces.add(array_api_compat.numpy)

    if single_namespace and len(namespaces) != 1:
        raise ValueError(
            f"Expected a single common namespace for array inputs, \
              but got: {[n.__name__ for n in namespaces]}"
        )

    (xp,) = namespaces

    return xp


def asarray(array, dtype=None, order=None, copy=None, *, xp=None):
    """Drop-in replacement for `np.asarray`.

    Memory layout parameter `order` is not exposed in the Array API standard.
    `order` is only enforced if the input array implementation
    is NumPy based, otherwise `order` is just silently ignored.

    The purpose of this helper is to make it possible to share code for data
    container validation without memory copies for both downstream use cases
    """
    if xp is None:
        xp = namespace_from_arrays(array)
    if xp.__name__ in {"numpy", "array_api_compat.numpy", "numpy.array_api"}:
        # Use NumPy API to support order
        if copy is True:
            array = np.array(array, order=order, dtype=dtype)
        else:
            array = np.asarray(array, order=order, dtype=dtype)

        # At this point array is a NumPy ndarray. We convert it to an array
        # container that is consistent with the input's namespace.
        return xp.asarray(array)
    else:
        return xp.asarray(array, dtype=dtype, copy=copy)

def asarray_namespace(*arrays):
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

    Examples
    --------
    >>> import numpy as np
    >>> x, y, xp = asarray_namespace([0, 1, 2], np.arange(3))
    >>> xp.__name__
    'array_api_compat.numpy'
    >>> x, y
    (array([0, 1, 2]]), array([0, 1, 2]))

    """
    arrays = list(arrays)  # probably not good
    xp = namespace_from_arrays(*arrays)

    for i, array in enumerate(arrays):
        arrays[i] = asarray(array, xp=xp)

    return *arrays, xp
