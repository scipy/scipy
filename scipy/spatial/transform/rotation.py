# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.

from __future__ import annotations


import numpy as np

from scipy._lib.deprecation import _sub_module_deprecation
import scipy.spatial.transform._rotation as cython_backend
import scipy.spatial.transform._rotation_array_api as array_api_backend
from scipy._lib._array_api import array_namespace, Array, is_numpy

__all__ = ["Rotation", "Slerp"]  # noqa: F822

# Fast path for numpy arrays: If quat is a numpy object, we call the Cython backend
backend_registry = {array_namespace(np.empty(0)): cython_backend}


class Rotation:
    def __init__(self, quat, normalize=True, copy=True, scalar_first=False):
        quat = self._sanitize_array(quat)
        # Legacy behavior for cython backend: Differentiate between single quat and batched quats
        xp = array_namespace(quat)
        self._single = quat.ndim == 1
        quat = xp.atleast_2d(quat)
        self._backend = backend_registry.get(xp, array_api_backend)
        self._quat: Array = self._backend.from_quat(
            quat, normalize=normalize, copy=copy, scalar_first=scalar_first
        )

    @classmethod
    def from_quat(cls, quat: Array, *, scalar_first=False) -> Rotation:
        return cls(quat, normalize=True, scalar_first=scalar_first)

    def as_quat(self, canonical=False, *, scalar_first=False):
        return self._backend.as_quat(
            self._quat, normalize=False, canonical=canonical, scalar_first=scalar_first
        )

    def as_matrix(self) -> Array:
        return self._backend.as_matrix(self._quat)

    def apply(self, points: Array, inverse=False) -> Array:
        return self._backend.apply(self._quat, points, inverse=inverse)

    def __getstate__(self):
        return self._xp.asarray(self._quat, dtype=float)

    def __setstate__(self, state):
        xp = array_namespace(state)
        self._backend = backend_registry.get(xp, array_api_backend)
        self._quat = xp.asarray(state, copy=True)

    @property
    def single(self):
        """Whether this instance represents a single rotation."""
        # TODO: Remove this once we properly support broadcasting with arbitrary
        # number of rotations
        return self._single

    def __bool__(self):
        """Comply with Python convention for objects to be True.

        Required because `Rotation.__len__()` is defined and not always truthy.
        """
        return True

    def __len__(self):
        """Number of rotations contained in this object.

        Multiple rotations can be stored in a single instance.

        Returns
        -------
        length : int
            Number of rotations stored in object.

        Raises
        ------
        TypeError if the instance was created as a single rotation.
        """
        if self._single:
            raise TypeError("Single rotation has no len().")

        return self._quat.shape[:-1].prod()

    def _sanitize_array(self, quat: Array) -> Array:
        xp = array_namespace(quat)
        quat = xp.asarray(quat)
        # Use float64 if available, else float32 (e.g. for jax)
        if quat.dtype == xp.float32 and not is_numpy(xp):
            dtype = xp.float32
        else:
            dtype = xp.result_type(xp.float32, xp.float64)
        # Always promote to float64 to make it compatible with the cython backend signatures
        return xp.asarray(quat, dtype=dtype)


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


# Register as pytree node for JAX if available to make Rotation compatible as input argument and
# return type for jit-compiled functions
try:
    from jax.tree_util import register_pytree_node

    def rot_unflatten(_, c):
        # Optimization: We do not want to call __init__ here because it would perform normalizations
        # twice. More importantly, it would call the non-jitted Array API backend and therefore
        # incur a significant performance hit
        r = Rotation.__new__(Rotation)
        r._backend = array_api_backend
        r._quat = c[0]
        # TODO: Remove this once we properly support broadcasting
        r._single = r._quat.ndim == 1
        return r

    register_pytree_node(Rotation, lambda v: ((v._quat,), None), rot_unflatten)
except ImportError:
    pass
