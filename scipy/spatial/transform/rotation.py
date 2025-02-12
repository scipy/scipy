# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.

from __future__ import annotations

from scipy._lib.deprecation import _sub_module_deprecation
from scipy._lib.array_api_compat import device
from ._rotation import Rotation as CythonRotation
from ._rotation_array_api import _as_quat, _as_matrix, _apply
from scipy._lib._array_api import array_namespace, is_numpy, Array, is_jax

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
        # Fast path for numpy arrays: If quat is a numpy object, we call the Cython version of the
        # class and forward all calls to the function's methods to super()
        self._use_cython = is_numpy(array_namespace(quat))
        if self._use_cython:
            super().__init__(
                quat, normalize=normalize, copy=copy, scalar_first=scalar_first
            )
            return
        self._quat = _as_quat(
            quat, normalize=normalize, copy=copy, scalar_first=scalar_first
        )

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
