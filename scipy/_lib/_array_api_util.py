"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""

import numpy as np

from scipy._lib import _pep440

_NUMPY_INCLUDES_ARRAY_API = _pep440.parse(np.__version__) >= _pep440.Version(
    "1.22.0"
)


def _get_namespace(*arrays):
    """
    Returns the module that implements Array API compatible functions or NumPy.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    *arrays : sequence of array_like
        Arrays to get the namespace from.

    Examples
    --------
    >>> import numpy as np
    >>> x = np.arange(6)
    >>> _get_namespace(x).__name__
    'numpy'

    >>> import numpy.array_api as xp
    >>> x = xp.arange(6)
    >>> _get_namespace(x).__name__
    'numpy.array_api'
    """
    # `arrays` contains one or more arrays
    # When there is no __array_namespace__,
    # we will use assume the default (NumPy) namespace.
    namespaces = {
        x.__array_namespace__() if hasattr(x, "__array_namespace__") else np
        for x in arrays
    }

    if len(namespaces) != 1:
        raise ValueError(
            f"Expected a single common namespace for array inputs, \
              but got: {[n.__name__ for n in namespaces]}"
        )

    (xp,) = namespaces

    return xp


def _concatenate(arrays, axis):
    """NumPy and Array API compatibility function for concatenating arrays."""
    xp = _get_namespace(*arrays)
    if xp is np:
        return xp.concatenate(arrays, axis=axis)
    else:
        return xp.concat(arrays, axis=axis)
