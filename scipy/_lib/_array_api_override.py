"""
Override functions from array_api_compat, for use by array-api-extra
and internally.

See also _array_api_compat_vendor.py
"""
import os

from collections.abc import Iterable, Iterator
from types import ModuleType
from typing import Any, TypeAlias

import numpy as np
import numpy.typing as npt

from scipy._lib import array_api_compat
import scipy._lib.array_api_compat.numpy as np_compat
from scipy._lib.array_api_compat import is_array_api_obj
from scipy._lib._sparse import issparse


Array: TypeAlias = Any  # To be changed to a Protocol later (see array-api#589)
ArrayLike: TypeAlias = Array | npt.ArrayLike

# To enable array API and strict array-like input validation
SCIPY_ARRAY_API: str | bool = os.environ.get("SCIPY_ARRAY_API", False)
# To control the default device - for use in the test suite only
SCIPY_DEVICE = os.environ.get("SCIPY_DEVICE", "cpu")


def _compliance_scipy(arrays: Iterable[ArrayLike]) -> Iterator[Array]:
    """Raise exceptions on known-bad subclasses. Discard 0-dimensional ArrayLikes
    and convert 1+-dimensional ArrayLikes to numpy.

    The following subclasses are not supported and raise and error:
    - `numpy.ma.MaskedArray`
    - `numpy.matrix`
    - NumPy arrays which do not have a boolean or numerical dtype
    - Any array-like which is neither array API compatible nor coercible by NumPy
    - Any array-like which is coerced by NumPy to an unsupported dtype
    """
    for array in arrays:
        if array is None:
            continue

        # this comes from `_util._asarray_validated`
        if issparse(array):
            msg = ('Sparse arrays/matrices are not supported by this function. '
                   'Perhaps one of the `scipy.sparse.linalg` functions '
                   'would work instead.')
            raise ValueError(msg)

        if isinstance(array, np.ma.MaskedArray):
            raise TypeError("Inputs of type `numpy.ma.MaskedArray` are not supported.")

        if isinstance(array, np.matrix):
            raise TypeError("Inputs of type `numpy.matrix` are not supported.")

        if isinstance(array, np.ndarray | np.generic):
            dtype = array.dtype
            if not (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bool_)):
                raise TypeError(f"An argument has dtype `{dtype!r}`; "
                                f"only boolean and numerical dtypes are supported.")

        if is_array_api_obj(array):
            yield array
        else:
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
            # Ignore 0-dimensional arrays, coherently with array-api-compat.
            # Raise if there are 1+-dimensional array-likes mixed with non-numpy
            # Array API objects.
            if array.ndim:
                yield array


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

    1. Check for the global switch: SCIPY_ARRAY_API.
    2. `_compliance_scipy` raise exceptions on known-bad subclasses. See
       its definition for more details.

    When the global switch is False, it defaults to the `numpy` namespace.
    In that case, there is no compliance check. This is a convenience to
    ease the adoption. Otherwise, arrays must comply with the new rules.
    """
    if not SCIPY_ARRAY_API:
        # here we could wrap the namespace if needed
        return np_compat

    api_arrays = list(_compliance_scipy(arrays))
    # In case of a mix of array API compliant arrays and scalars, return
    # the array API namespace. If there are only ArrayLikes (e.g. lists),
    # return NumPy (wrapped by array-api-compat).
    if api_arrays:
        return array_api_compat.array_namespace(*api_arrays)
    return np_compat
