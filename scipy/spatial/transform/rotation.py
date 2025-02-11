# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.

from __future__ import annotations

from scipy._lib.deprecation import _sub_module_deprecation
from scipy._lib.array_api_compat import device
from ._rotation import Rotation as CythonRotation
from ._rotation_array_api import _as_quat, _as_matrix, _apply
from scipy._lib._array_api import array_namespace, is_numpy, Array

try:
    import jax  # Must succeed if we get passed a jax.numpy array
except ImportError:
    jax = None

__all__ = [  # noqa: F822
    "Rotation",
    "Slerp",
]

from functools import wraps
from ._rotation import Rotation as CythonRotation


def _maybe_cython(f):
    cython_fn = getattr(CythonRotation, f.__name__)

    @wraps(f)
    def wrapper(self, *args, **kwargs):
        if self._use_cython:
            return cython_fn(self, *args, **kwargs)
        return f(self, *args, **kwargs)

    return wrapper


class Rotation(CythonRotation):
    def __init__(self, quat, normalize=True, copy=True, scalar_first=False):
        xp = array_namespace(quat)
        # Fast path for numpy arrays: If quat is a numpy object, we call the Cython version of the
        # class and forward all calls to the function's methods to super()
        self._use_cython = is_numpy(xp)
        self._xp = xp
        if self._use_cython:
            super().__init__(
                quat, normalize=normalize, copy=copy, scalar_first=scalar_first
            )
            return
        dtype = xp.float64 if quat.dtype not in (xp.float32, xp.float64) else quat.dtype
        quat = xp.asarray(quat, dtype=dtype)
        if quat.shape[-1] != 4 or quat.shape[0] == 0:
            raise ValueError(
                f"Expected `quat` to have shape (..., 4), got {quat.shape}."
            )
        scalar_first = xp.asarray(scalar_first)
        normalize = xp.asarray(normalize)
        copy = xp.asarray(copy)
        # All operations are done non-branching to enable JIT-compilation
        quat = xp.where(scalar_first, xp.roll(quat, -1, axis=-1), quat)
        quat = xp.where(normalize | copy, xp.asarray(quat, copy=True), quat)
        quat_norm = xp.linalg.vector_norm(quat, axis=-1, keepdims=True)

        # TODO: This is a deviation from the behavior of the cython version because JIT code needs
        # to be strictly non-branching. If we have zero quaternions, we return nan values
        # if xp.any(quat_norm == 0):
        #     raise ValueError("Found zero norm quaternions in `quat`.")
        quat_norm = xp.where(
            quat_norm == 0, xp.asarray([xp.nan], device=device(quat)), quat_norm
        )
        self._quat = xp.where(normalize, quat / quat_norm, quat)

    @classmethod
    def from_quat(cls, quat: Array, *, scalar_first=False) -> Rotation:
        return cls(quat, normalize=True, scalar_first=scalar_first)

    @_maybe_cython
    def as_quat(self, canonical=False, *, scalar_first=False):
        return _as_quat(self._quat, canonical=canonical, scalar_first=scalar_first)

    @_maybe_cython
    def as_matrix(self) -> Array:
        return _as_matrix(self._quat)

    @_maybe_cython
    def apply(self, points: Array) -> Array:
        return _apply(self._quat, points)

    @_maybe_cython
    def __getstate__(self):
        return self._xp.asarray(self._quat, dtype=float), None

    def __setstate__(self, state):
        quat, single = state
        if single is not None:  # Cython version
            self._use_cython = True
            self._xp = array_namespace(quat)
            return super().__setstate__(state)
        self._xp = array_namespace(quat)
        self._quat = self._xp.asarray(quat, copy=True)


def __dir__():
    return __all__


def __getattr__(name):
    return _sub_module_deprecation(
        sub_package="spatial.transform",
        module="rotation",
        private_modules=["_rotation"],
        all=__all__,
        attribute=name,
    )
