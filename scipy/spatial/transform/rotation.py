# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.

from __future__ import annotations


import numpy as np

from scipy._lib.deprecation import _sub_module_deprecation
import scipy.spatial.transform._cython_backend as cython_backend
import scipy.spatial.transform._array_api_backend as array_api_backend
from scipy._lib._array_api import array_namespace, Array, is_numpy, ArrayLike
import scipy._lib.array_api_extra as xpx


__all__ = ["Rotation", "Slerp"]  # noqa: F822

# Fast path for numpy arrays: If quat is a numpy object, we call the Cython backend
backend_registry = {array_namespace(np.empty(0)): cython_backend}


class Rotation:
    def __init__(
        self,
        quat: ArrayLike,
        normalize: bool = True,
        copy: bool = True,
        scalar_first: bool = False,
    ):
        quat = self._to_array(quat)
        xp = array_namespace(quat)
        # Legacy behavior for cython backend: Differentiate between single quat and batched quats
        # We only use this for the cython backend. The Array API backend uses broadcasting by
        # default and hence returns the correct shape without additional logic
        self._single = quat.ndim == 1 and is_numpy(xp)
        if self._single:
            quat = xpx.atleast_nd(quat, ndim=2, xp=xp)
        self._backend = backend_registry.get(xp, array_api_backend)
        self._quat: Array = self._backend.from_quat(
            quat, normalize=normalize, copy=copy, scalar_first=scalar_first
        )

    @classmethod
    def from_quat(cls, quat: ArrayLike, *, scalar_first: bool = False) -> Rotation:
        return cls(quat, normalize=True, scalar_first=scalar_first)

    @classmethod
    def from_euler(
        cls, seq: str, angles: ArrayLike, *, degrees: bool = False
    ) -> Rotation:
        backend = backend_registry.get(array_namespace(angles), array_api_backend)
        quat = backend.from_euler(seq, angles, degrees=degrees)
        return cls(quat, normalize=False, copy=False)

    @classmethod
    def from_matrix(cls, matrix: ArrayLike) -> Rotation:
        backend = backend_registry.get(array_namespace(matrix), array_api_backend)
        quat = backend.from_matrix(matrix)
        return cls(quat, normalize=False, copy=False)

    def as_quat(self, canonical=False, *, scalar_first=False):
        quat = self._backend.as_quat(
            self._quat, canonical=canonical, scalar_first=scalar_first
        )
        if self._single:
            return quat[0, ...]
        return quat

    def as_matrix(self) -> Array:
        matrix = self._backend.as_matrix(self._quat)
        if self._single:
            return matrix[0, ...]
        return matrix

    def apply(self, points: Array, inverse: bool = False) -> Array:
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

    def inv(self):
        q_inv = self._backend.inv(self._quat)
        if self._single:
            q_inv = q_inv[0]
        return Rotation(q_inv, normalize=False, copy=False)

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
        # TODO: Should we replace this with a shape instead?
        if self._single:
            raise TypeError("Single rotation has no len().")

        return self._quat.shape[:-1].prod()

    def __mul__(self, other):
        if not _broadcastable(self._quat.shape, other._quat.shape):
            raise ValueError(
                "Expected equal number of rotations in both or a single "
                f"rotation in either object, got {self._quat.shape[:-1]} rotations in "
                f"first and {other._quat.shape[:-1]} rotations in second object."
            )
        q_composed = self._backend.compose_quat(self._quat, other._quat)
        if self._single:
            q_composed = q_composed[0]
        return Rotation(q_composed, normalize=True, copy=False)

    def _to_array(self, quat: ArrayLike) -> Array:
        """Convert the quaternion to an array.

        The return array dtype follows the following rules:
        - If quat is an ArrayLike or NumPy array, we always promote to float64
        - If quat is an Array from frameworks other than NumPy, we preserve the dtype if it is
          float32. Otherwise, we promote to the result type of combining float32 and float64

        The first rule is required by the cython backend signatures that expect cython.double views.
        The second rule is necessary to promote non-floating arrays to the correct type in
        frameworks that may not support double precision (e.g. jax by default).
        """
        xp = array_namespace(quat)
        quat = xp.asarray(quat)
        # TODO: Remove this once we properly support broadcasting
        if quat.ndim not in (1, 2) or quat.shape[-1] != 4 or quat.shape[0] == 0:
            raise ValueError(f"Expected `quat` to have shape (N, 4), got {quat.shape}.")

        if quat.dtype == xp.float32 and not is_numpy(xp):
            dtype = xp.float32
        else:
            dtype = xp.result_type(xp.float32, xp.float64)
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


def _broadcastable(shape_a: tuple[int, ...], shape_b: tuple[int, ...]) -> bool:
    """Check if two shapes are broadcastable."""
    return all(
        (m == n) or (m == 1) or (n == 1) for m, n in zip(shape_a[::-1], shape_b[::-1])
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
