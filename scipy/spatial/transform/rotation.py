# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.

from scipy._lib.deprecation import _sub_module_deprecation
from ._rotation import Rotation as CRotation
from ._rotation_array_api import _as_quat, _as_matrix
from scipy._lib._array_api import array_namespace, is_numpy, is_jax
import time

try:
    import jax  # Must succeed if we get passed a jax.numpy array
except ImportError:
    jax = None

__all__ = [  # noqa: F822
    "Rotation",
    "Slerp",
]


class Rotation(CRotation):
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
        # TODO: Do we want this dtype check?
        dtype = (
            xp.float64 if not is_jax(xp) or jax.config.jax_enable_x64 else xp.float32
        )
        quat = xp.asarray(quat, dtype=dtype)
        if quat.shape[-1] != 4 or quat.shape[0] == 0:
            raise ValueError(
                f"Expected `quat` to have shape (..., 4), got {quat.shape}."
            )
        scalar_first = xp.asarray(scalar_first, device=quat.device)
        normalize = xp.asarray(normalize, device=quat.device)
        copy = xp.asarray(copy, device=quat.device)
        # All operations are done non-branching to enable JIT-compilation
        quat = xp.where(scalar_first, xp.roll(quat, -1, axis=-1), quat)
        quat = xp.where(normalize | copy, xp.asarray(quat, copy=True), quat)
        quat_norm = xp.linalg.vector_norm(quat, axis=-1, keepdims=True)
        if xp.any(quat_norm == 0):
            raise ValueError("Found zero norm quaternions in `quat`.")
        self._quat = xp.where(normalize, quat / quat_norm, quat)

    @classmethod
    def from_quat(cls, quat, *, scalar_first=False):
        return cls(quat, normalize=True, scalar_first=scalar_first)

    def as_quat(self, canonical=False, *, scalar_first=False):
        if self._use_cython:
            return super().as_quat(canonical, scalar_first=scalar_first)
        return _as_quat(self._quat, canonical=canonical, scalar_first=scalar_first)

    def as_matrix(self):
        if self._use_cython:
            return super().as_matrix()
        return _as_matrix(self._quat)

    def __getstate__(self):
        if self._use_cython:
            return super().__getstate__()
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
