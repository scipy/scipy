"""Array-API backend for scipy.spatial.distance."""

from types import ModuleType
from typing import Any

import numpy as np

from scipy._lib._array_api import (
    Array,
    ArrayLike,
    _asarray,
    array_namespace,
    is_lazy_array,
    is_numpy,
    xp_promote,
)

# Functions in the backend should be batched by default, operating over the last axis
# and broadcasting over arbitrary leading dimensions. This allows us to keep one
# implementation for the single-pair and pdist/cdist cases


def _validate_vector(u: ArrayLike, dtype: Any | None = None) -> Array:
    try:
        u = _asarray(u, dtype=dtype)
    except TypeError:
        # String/object arrays (e.g. hamming on byte strings) have no array-API
        # namespace. Fall back to NumPy, which handles them.
        u = np.asarray(u, dtype=dtype)
    if u.ndim == 1:
        return u
    raise ValueError("Input vector should be 1-D.")


def _validate_weights(w: ArrayLike, xp: ModuleType | None = None) -> Array:
    """Reject negatives along the last axis.

    Note:
        If w is a numpy array, we also cast to float64. If w is lazy, we fill the
        weights with NaN.
    """
    xp = array_namespace(w) if xp is None else xp
    w = _promote(w, xp=xp)
    invalid = w < 0
    if is_lazy_array(w):
        any_invalid = xp.any(invalid, axis=-1, keepdims=True)
        return xp.where(any_invalid, xp.nan, w)
    if xp.any(invalid):
        raise ValueError("Input weights should be all non-negative")
    return w


def _promote(x: ArrayLike, xp: ModuleType) -> Array:
    """Promote array to float64 for numpy, else according to the Array API spec.

    The return array dtype follows the following rules:
    - If x is an ArrayLike or NumPy array, we always promote to float64
    - If x is an Array from frameworks other than NumPy, we preserve the precision
      of the input array dtype.
    """
    if is_numpy(xp):
        return _asarray(x, order="C", xp=xp, dtype=xp.float64)
    return xp_promote(x, force_floating=True, xp=xp)


def minkowski(u: ArrayLike, v: ArrayLike, p: float = 2, w: ArrayLike | None = None):
    xp = array_namespace(u, v)
    u = _promote(u, xp=xp)
    v = _promote(v, xp=xp)
    if p <= 0:
        raise ValueError("p must be greater than 0")
    u_v = u - v
    if w is not None:
        w = _validate_weights(w, xp=xp)
        if p == 1:
            root_w = w
        elif p == 2:
            root_w = xp.sqrt(w)  # better precision and speed
        elif p == np.inf:
            root_w = xp.astype(w != 0, w.dtype)
        else:
            root_w = w ** (1 / p)
        u_v = root_w * u_v
    return xp.linalg.vector_norm(u_v, ord=p, axis=-1)
