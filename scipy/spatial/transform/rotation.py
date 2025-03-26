# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.

from __future__ import annotations

from typing import Iterable

import numpy as np

from scipy._lib.deprecation import _sub_module_deprecation
import scipy.spatial.transform._cython_backend as cython_backend
import scipy.spatial.transform._array_api_backend as array_api_backend
from scipy._lib._array_api import array_namespace, Array, is_numpy, ArrayLike, is_jax
from scipy._lib.array_api_compat import device
import scipy._lib.array_api_extra as xpx
from scipy._lib._util import _transition_to_rng


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
    def from_matrix(cls, matrix: ArrayLike) -> Rotation:
        backend = backend_registry.get(array_namespace(matrix), array_api_backend)
        quat = backend.from_matrix(matrix)
        return cls(quat, normalize=False, copy=False)

    @classmethod
    def from_rotvec(cls, rotvec: ArrayLike, degrees: bool = False) -> Rotation:
        backend = backend_registry.get(array_namespace(rotvec), array_api_backend)
        quat = backend.from_rotvec(rotvec, degrees=degrees)
        return cls(quat, normalize=False, copy=False)

    @classmethod
    def from_mrp(cls, mrp: ArrayLike) -> Rotation:
        backend = backend_registry.get(array_namespace(mrp), array_api_backend)
        quat = backend.from_mrp(mrp)
        return cls(quat, normalize=False, copy=False)

    @classmethod
    def from_euler(cls, seq: str, angles: ArrayLike, degrees: bool = False) -> Rotation:
        backend = backend_registry.get(array_namespace(angles), array_api_backend)
        quat = backend.from_euler(seq, angles, degrees=degrees)
        return cls(quat, normalize=False, copy=False)

    @classmethod
    def from_davenport(
        cls,
        axes: ArrayLike,
        order: str,
        angles: ArrayLike | float,
        degrees: bool = False,
    ) -> Rotation:
        backend = backend_registry.get(array_namespace(axes), array_api_backend)
        quat = backend.from_davenport(axes, order, angles, degrees)
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

    def as_rotvec(self, degrees: bool = False) -> Array:
        rotvec = self._backend.as_rotvec(self._quat, degrees=degrees)
        if self._single:
            return rotvec[0, ...]
        return rotvec

    def as_mrp(self) -> Array:
        mrp = self._backend.as_mrp(self._quat)
        if self._single:
            return mrp[0, ...]
        return mrp

    def as_euler(self, seq: str, degrees: bool = False) -> Array:
        euler = self._backend.as_euler(self._quat, seq, degrees=degrees)
        if self._single:
            return euler[0, ...]
        return euler

    def as_davenport(self, axes: ArrayLike, order: str, degrees: bool = False) -> Array:
        """Represent as Davenport angles.

        Any orientation can be expressed as a composition of 3 elementary
        rotations.

        For both Euler angles and Davenport angles, consecutive axes must
        be are orthogonal (``axis2`` is orthogonal to both ``axis1`` and
        ``axis3``). For Euler angles, there is an additional relationship
        between ``axis1`` or ``axis3``, with two possibilities:

            - ``axis1`` and ``axis3`` are also orthogonal (asymmetric sequence)
            - ``axis1 == axis3`` (symmetric sequence)

        For Davenport angles, this last relationship is relaxed [1]_, and only
        the consecutive orthogonal axes requirement is maintained.

        A slightly modified version of the algorithm from [2]_ has been used to
        calculate Davenport angles for the rotation about a given sequence of
        axes.

        Davenport angles, just like Euler angles, suffer from the problem of
        gimbal lock [3]_, where the representation loses a degree of freedom
        and it is not possible to determine the first and third angles
        uniquely. In this case, a warning is raised, and the third angle is set
        to zero. Note however that the returned angles still represent the
        correct rotation.

        Parameters
        ----------
        axes : array_like, shape (3,) or ([1 or 2 or 3], 3)
            Axis of rotation, if one dimensional. If two dimensional, describes the
            sequence of axes for rotations, where each axes[i, :] is the ith
            axis. If more than one axis is given, then the second axis must be
            orthogonal to both the first and third axes.
        order : string
            If it belongs to the set {'e', 'extrinsic'}, the sequence will be
            extrinsic. If if belongs to the set {'i', 'intrinsic'}, sequence
            will be treated as intrinsic.
        degrees : boolean, optional
            Returned angles are in degrees if this flag is True, else they are
            in radians. Default is False.

        Returns
        -------
        angles : ndarray, shape (3,) or (N, 3)
            Shape depends on shape of inputs used to initialize object.
            The returned angles are in the range:

            - First angle belongs to [-180, 180] degrees (both inclusive)
            - Third angle belongs to [-180, 180] degrees (both inclusive)
            - Second angle belongs to a set of size 180 degrees,
                given by: ``[-abs(lambda), 180 - abs(lambda)]``, where ``lambda``
                is the angle between the first and third axes.

        References
        ----------
        .. [1] Shuster, Malcolm & Markley, Landis. (2003). Generalization of
                the Euler Angles. Journal of the Astronautical Sciences. 51. 123-132. 10.1007/BF03546304.
        .. [2] Bernardes E, Viollet S (2022) Quaternion to Euler angles
                conversion: A direct, general and computationally efficient method.
                PLoS ONE 17(11): e0276302. 10.1371/journal.pone.0276302
        .. [3] https://en.wikipedia.org/wiki/Gimbal_lock#In_applied_mathematics

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Davenport angles are a generalization of Euler angles, when we use the
        canonical basis axes:

        >>> ex = [1, 0, 0]
        >>> ey = [0, 1, 0]
        >>> ez = [0, 0, 1]

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True)
        array([90.,  0.,  0.])
        >>> r.as_euler('zxy', degrees=True)
        array([90.,  0.,  0.])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True).shape
        (3,)

        Represent a stack of single rotation:

        >>> r = R.from_rotvec([[0, 0, np.pi/2]])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True)
        array([[90.,  0.,  0.]])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True).shape
        (1, 3)

        Represent multiple rotations in a single object:

        >>> r = R.from_rotvec([
        ... [0, 0, 90],
        ... [45, 0, 0]], degrees=True)
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True)
        array([[90.,  0.,  0.],
                [ 0., 45.,  0.]])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True).shape
        (2, 3)
        """
        davenport = self._backend.as_davenport(self._quat, axes, order, degrees)
        if self._single:
            return davenport[0, ...]
        return davenport

    @classmethod
    @_transition_to_rng("random_state", position_num=2)
    def random(cls, num: int | None = None, rng: np.random.Generator | None = None):
        # DECISION: How do we handle random numbers in other frameworks?
        # TODO: The array API does not have a unified random interface. This method only creates
        # numpy arrays. If we do want to support other frameworks, we need a way to handle other rng
        # implementations.
        sample = cython_backend.random(num, rng)
        return cls(sample, normalize=True, copy=False)

    @classmethod
    def identity(cls, num: int | None = None) -> Rotation:
        return cls(cython_backend.identity(num), normalize=False, copy=False)

    def inv(self):
        """Invert this rotation.

        Composition of a rotation with its inverse results in an identity
        transformation.

        Returns
        -------
        inverse : `Rotation` instance
            Object containing inverse of the rotations in the current instance.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Inverting a single rotation:

        >>> p = R.from_euler('z', 45, degrees=True)
        >>> q = p.inv()
        >>> q.as_euler('zyx', degrees=True)
        array([-45.,   0.,   0.])

        Inverting multiple rotations:

        >>> p = R.from_rotvec([[0, 0, np.pi/3], [-np.pi/4, 0, 0]])
        >>> q = p.inv()
        >>> q.as_rotvec()
        array([[-0.        , -0.        , -1.04719755],
               [ 0.78539816, -0.        , -0.        ]])

        """
        q_inv = self._backend.inv(self._quat)
        if self._single:
            q_inv = q_inv[0, ...]
        return Rotation(q_inv, normalize=False, copy=False)

    def magnitude(self):
        magnitude = self._backend.magnitude(self._quat)
        if self._single:
            return magnitude[0, ...]
        return magnitude

    def approx_equal(
        self, other: Rotation, atol: float | None = None, degrees: bool = False
    ):
        return self._backend.approx_equal(
            self._quat, other._quat, atol=atol, degrees=degrees
        )

    def mean(self, weights: ArrayLike | None = None) -> Rotation:
        return Rotation(
            self._backend.mean(self._quat, weights=weights), normalize=False
        )

    def reduce(
        self,
        left: Rotation | None = None,
        right: Rotation | None = None,
        return_indices: bool = False,
    ) -> Rotation:
        left = left.as_quat() if left is not None else None
        right = right.as_quat() if right is not None else None
        reduced, left_idx, right_idx = self._backend.reduce(
            self._quat, left=left, right=right
        )
        rot = Rotation(reduced, normalize=False, copy=False)
        if return_indices:
            left_idx = left_idx if left is not None else None
            right_idx = right_idx if right is not None else None
            return rot, left_idx, right_idx
        return rot

    def apply(self, points: Array, inverse: bool = False) -> Array:
        points = array_namespace(self._quat).asarray(
            points,
            device=device(self._quat),
            dtype=array_api_backend.atleast_f32(self._quat),
        )
        result = self._backend.apply(self._quat, points, inverse=inverse)
        if self._single and points.ndim == 1:
            return result[0, ...]
        return result

    @classmethod
    def align_vectors(
        cls,
        a: Array,
        b: Array,
        weights: Array | None = None,
        return_sensitivity: bool = False,
    ) -> tuple[Rotation, float] | tuple[Rotation, float, Array]:
        backend = backend_registry.get(array_namespace(a), array_api_backend)
        q, rssd, sensitivity = backend.align_vectors(a, b, weights, return_sensitivity)
        if return_sensitivity:
            return cls(q, normalize=False, copy=False), rssd, sensitivity
        return cls(q, normalize=False, copy=False), rssd

    @classmethod
    def concatenate(cls, rotations: Rotation | Iterable[Rotation]) -> Rotation:
        """Concatenate a sequence of `Rotation` objects into a single object.

        This is useful if you want to, for example, take the mean of a set of
        rotations and need to pack them into a single object to do so.

        Parameters
        ----------
        rotations : sequence of `Rotation` objects
            The rotations to concatenate. If a single `Rotation` object is
            passed in, a copy is returned.

        Returns
        -------
        concatenated : `Rotation` instance
            The concatenated rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> r1 = R.from_rotvec([0, 0, 1])
        >>> r2 = R.from_rotvec([0, 0, 2])
        >>> rc = R.concatenate([r1, r2])
        >>> rc.as_rotvec()
        array([[0., 0., 1.],
                [0., 0., 2.]])
        >>> rc.mean().as_rotvec()
        array([0., 0., 1.5])

        Concatenation of a split rotation recovers the original object.

        >>> rs = [r for r in rc]
        >>> R.concatenate(rs).as_rotvec()
        array([[0., 0., 1.],
                [0., 0., 2.]])

        Note that it may be simpler to create the desired rotations by passing
        in a single list of the data during initialization, rather then by
        concatenating:

        >>> R.from_rotvec([[0, 0, 1], [0, 0, 2]]).as_rotvec()
        array([[0., 0., 1.],
                [0., 0., 2.]])

        Notes
        -----
        .. versionadded:: 1.8.0
        """
        if isinstance(rotations, Rotation):
            return cls(rotations.as_quat(), normalize=False, copy=True)
        if not all(isinstance(x, Rotation) for x in rotations):
            raise TypeError("input must contain Rotation objects only")

        xp = array_namespace(rotations[0].as_quat())
        quats = xp.concat(
            [xpx.atleast_nd(x.as_quat(), ndim=2, xp=xp) for x in rotations]
        )
        return cls(quats, normalize=False)

    def __getitem__(self, indexer) -> Rotation:
        if self._single or self._quat.ndim == 1:
            raise TypeError("Single rotation is not subscriptable.")
        is_array = isinstance(indexer, type(self._quat))
        # Masking is only specified in the Array API when the array is the sole index
        # TODO: Make access to xp more efficient
        # TODO: This special case handling is mainly a result of Array API limitations. Ideally we
        # would get rid of them altogether and converge to [indexer, ...] indexing.
        xp = array_namespace(self._quat)
        if is_array and indexer.dtype == xp.bool:
            return Rotation(self._quat[indexer], normalize=False)
        if is_array and (indexer.dtype == xp.int64 or indexer.dtype == xp.int32):
            # Array API limitation: Integer index arrays are only allowed with integer indices
            all_ind = xp.arange(4)
            indexer = xp.reshape(indexer, (indexer.shape[0], 1))
            return Rotation(self._quat[indexer, all_ind], normalize=False)
        return Rotation(self._quat[indexer, ...], normalize=False)

    def __setitem__(self, indexer, value: Rotation):
        if self._single or self._quat.ndim == 1:
            raise TypeError("Single rotation is not subscriptable.")

        if not isinstance(value, Rotation):
            raise TypeError("value must be a Rotation object")

        self._quat = self._backend.setitem(self._quat, value.as_quat(), indexer)

    def __getstate__(self):
        return (self._quat, self._single)

    def __setstate__(self, state):
        quat, single = state
        xp = array_namespace(quat)
        self._backend = backend_registry.get(xp, array_api_backend)
        self._quat = xp.asarray(quat, copy=True)
        self._single = single

    @property
    def single(self):
        """Whether this instance represents a single rotation."""
        # TODO: Remove this once we properly support broadcasting with arbitrary
        # number of rotations
        return self._single or self._quat.ndim == 1

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
        if self._single or self._quat.ndim == 1:
            raise TypeError("Single rotation has no len().")
        return self._quat.shape[0]

    def __mul__(self, other):
        """Compose this rotation with the other.

        If `p` and `q` are two rotations, then the composition of 'q followed
        by p' is equivalent to `p * q`. In terms of rotation matrices,
        the composition can be expressed as
        ``p.as_matrix() @ q.as_matrix()``.

        Parameters
        ----------
        other : `Rotation` instance
            Object containing the rotations to be composed with this one. Note
            that rotation compositions are not commutative, so ``p * q`` is
            generally different from ``q * p``.

        Returns
        -------
        composition : `Rotation` instance
            This function supports composition of multiple rotations at a time.
            The following cases are possible:

            - Either ``p`` or ``q`` contains a single rotation. In this case
              `composition` contains the result of composing each rotation in
              the other object with the single rotation.
            - Both ``p`` and ``q`` contain ``N`` rotations. In this case each
              rotation ``p[i]`` is composed with the corresponding rotation
              ``q[i]`` and `output` contains ``N`` rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Composition of two single rotations:

        >>> p = R.from_quat([0, 0, 1, 1])
        >>> q = R.from_quat([1, 0, 0, 1])
        >>> p.as_matrix()
        array([[ 0., -1.,  0.],
               [ 1.,  0.,  0.],
               [ 0.,  0.,  1.]])
        >>> q.as_matrix()
        array([[ 1.,  0.,  0.],
               [ 0.,  0., -1.],
               [ 0.,  1.,  0.]])
        >>> r = p * q
        >>> r.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])

        Composition of two objects containing equal number of rotations:

        >>> p = R.from_quat([[0, 0, 1, 1], [1, 0, 0, 1]])
        >>> q = R.from_rotvec([[np.pi/4, 0, 0], [-np.pi/4, 0, np.pi/4]])
        >>> p.as_quat()
        array([[0.        , 0.        , 0.70710678, 0.70710678],
               [0.70710678, 0.        , 0.        , 0.70710678]])
        >>> q.as_quat()
        array([[ 0.38268343,  0.        ,  0.        ,  0.92387953],
               [-0.37282173,  0.        ,  0.37282173,  0.84971049]])
        >>> r = p * q
        >>> r.as_quat()
        array([[ 0.27059805,  0.27059805,  0.65328148,  0.65328148],
               [ 0.33721128, -0.26362477,  0.26362477,  0.86446082]])

        """
        if not array_api_backend.broadcastable(self._quat.shape, other._quat.shape):
            raise ValueError(
                "Expected equal number of rotations in both or a single "
                f"rotation in either object, got {self._quat.shape[:-1]} rotations in "
                f"first and {other._quat.shape[:-1]} rotations in second object."
            )
        quat = self._backend.compose_quat(self._quat, other._quat)
        if self._single and other._single:
            quat = quat[0]
        return Rotation(quat, normalize=True, copy=False)

    def __pow__(self, n: float, modulus: None = None) -> Rotation:
        """Compose this rotation with itself `n` times.

        Composition of a rotation ``p`` with itself can be extended to
        non-integer ``n`` by considering the power ``n`` to be a scale factor
        applied to the angle of rotation about the rotation's fixed axis. The
        expression ``q = p ** n`` can also be expressed as
        ``q = Rotation.from_rotvec(n * p.as_rotvec())``.

        If ``n`` is negative, then the rotation is inverted before the power
        is applied. In other words, ``p ** -abs(n) == p.inv() ** abs(n)``.

        Parameters
        ----------
        n : float
            The number of times to compose the rotation with itself.
        modulus : None
            This overridden argument is not applicable to Rotations and must be
            ``None``.

        Returns
        -------
        power : `Rotation` instance
            If the input Rotation ``p`` contains ``N`` multiple rotations, then
            the output will contain ``N`` rotations where the ``i`` th rotation
            is equal to ``p[i] ** n``

        Notes
        -----
        For example, a power of 2 will double the angle of rotation, and a
        power of 0.5 will halve the angle. There are three notable cases: if
        ``n == 1`` then the original rotation is returned, if ``n == 0``
        then the identity rotation is returned, and if ``n == -1`` then
        ``p.inv()`` is returned.

        Note that fractional powers ``n`` which effectively take a root of
        rotation, do so using the shortest path smallest representation of that
        angle (the principal root). This means that powers of ``n`` and ``1/n``
        are not necessarily inverses of each other. For example, a 0.5 power of
        a +240 degree rotation will be calculated as the 0.5 power of a -120
        degree rotation, with the result being a rotation of -60 rather than
        +120 degrees.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Raising a rotation to a power:

        >>> p = R.from_rotvec([1, 0, 0])
        >>> q = p ** 2
        >>> q.as_rotvec()
        array([2., 0., 0.])
        >>> r = p ** 0.5
        >>> r.as_rotvec()
        array([0.5, 0., 0.])

        Inverse powers do not necessarily cancel out:

        >>> p = R.from_rotvec([0, 0, 120], degrees=True)
        >>> ((p ** 2) ** 0.5).as_rotvec(degrees=True)
        array([  -0.,   -0., -60.])

        """
        if modulus is not None:
            raise NotImplementedError("modulus not supported")
        quat = self._backend.pow(self._quat, n)
        if self._single:
            quat = quat[0]
        return Rotation(quat, normalize=False, copy=False)

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
        if quat.ndim not in (1, 2) or quat.shape[-1] != 4:
            raise ValueError(f"Expected `quat` to have shape (N, 4), got {quat.shape}.")

        # TODO: Do we always want to promote to float64 for NumPy? This is consistent with the old
        # implementation, but it might make more sense to preserve float32 if passed in by the user.
        # This would make the behavior more consistent with the Array API backend, but requires
        # changes in the cython backend.
        if is_numpy(xp):
            dtype = xp.float64
        else:
            dtype = array_api_backend.atleast_f32(quat)
        return xp.asarray(quat, dtype=dtype)

    def __repr__(self):
        m = f"{np.asarray(self.as_matrix())!r}".splitlines()
        # bump indent (+21 characters)
        m[1:] = [" " * 21 + m[i] for i in range(1, len(m))]
        return "Rotation.from_matrix(" + "\n".join(m) + ")"


class Slerp:
    """Spherical Linear Interpolation of Rotations.

    The interpolation between consecutive rotations is performed as a rotation
    around a fixed axis with a constant angular velocity [1]_. This ensures
    that the interpolated rotations follow the shortest path between initial
    and final orientations.

    Parameters
    ----------
    times : array_like, shape (N,)
        Times of the known rotations. At least 2 times must be specified.
    rotations : `Rotation` instance
        Rotations to perform the interpolation between. Must contain N
        rotations.

    Methods
    -------
    __call__

    See Also
    --------
    Rotation

    Notes
    -----
    .. versionadded:: 1.2.0

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Slerp#Quaternion_Slerp

    Examples
    --------
    >>> from scipy.spatial.transform import Rotation as R
    >>> from scipy.spatial.transform import Slerp

    Setup the fixed keyframe rotations and times:

    >>> key_rots = R.random(5, random_state=2342345)
    >>> key_times = [0, 1, 2, 3, 4]

    Create the interpolator object:

    >>> slerp = Slerp(key_times, key_rots)

    Interpolate the rotations at the given times:

    >>> times = [0, 0.5, 0.25, 1, 1.5, 2, 2.75, 3, 3.25, 3.60, 4]
    >>> interp_rots = slerp(times)

    The keyframe rotations expressed as Euler angles:

    >>> key_rots.as_euler('xyz', degrees=True)
    array([[ 14.31443779, -27.50095894,  -3.7275787 ],
           [ -1.79924227, -24.69421529, 164.57701743],
           [146.15020772,  43.22849451, -31.34891088],
           [ 46.39959442,  11.62126073, -45.99719267],
           [-88.94647804, -49.64400082, -65.80546984]])

    The interpolated rotations expressed as Euler angles. These agree with the
    keyframe rotations at both endpoints of the range of keyframe times.

    >>> interp_rots.as_euler('xyz', degrees=True)
    array([[  14.31443779,  -27.50095894,   -3.7275787 ],
           [   4.74588574,  -32.44683966,   81.25139984],
           [  10.71094749,  -31.56690154,   38.06896408],
           [  -1.79924227,  -24.69421529,  164.57701743],
           [  11.72796022,   51.64207311, -171.7374683 ],
           [ 146.15020772,   43.22849451,  -31.34891088],
           [  68.10921869,   20.67625074,  -48.74886034],
           [  46.39959442,   11.62126073,  -45.99719267],
           [  12.35552615,    4.21525086,  -64.89288124],
           [ -30.08117143,  -19.90769513,  -78.98121326],
           [ -88.94647804,  -49.64400082,  -65.80546984]])

    """

    def __init__(self, times: ArrayLike, rotations: Rotation):
        if not isinstance(rotations, Rotation):
            raise TypeError("`rotations` must be a `Rotation` instance.")
        if rotations.single or len(rotations) <= 1:
            raise ValueError("`rotations` must be a sequence of at least 2 rotations.")
        q = rotations.as_quat()
        xp = array_namespace(q)
        times = xp.asarray(times, device=device(q), dtype=q.dtype)
        if times.ndim != 1:
            raise ValueError(
                "Expected times to be specified in a 1 dimensional array, got "
                f"{times.ndim} dimensions."
            )
        if times.shape[0] != len(rotations):
            raise ValueError(
                "Expected number of rotations to be equal to number of "
                f"timestamps given, got {len(rotations)} rotations and "
                f"{times.shape[0]} timestamps."
            )
        self.times = times
        # TODO: Replace with xp.diff once we upgrade to Array API 2024.12
        self.timedelta = times[1:] - times[:-1]

        # We cannot check for values for jit compiled code, so we cannot raise an error on timedelta
        # < 0 in jax. Instead, we set timedelta to nans
        # DECISION: Do we want to introduce this special case for jax or should all implementations
        # set the timedelta to nans?
        if is_jax(xp):
            mask = xp.any(self.timedelta <= 0)
            self.timedelta = xp.where(mask, xp.asarray(xp.nan), self.timedelta)
            self.times = xp.where(mask, xp.asarray(xp.nan), self.times)
        elif xp.any(self.timedelta <= 0):
            raise ValueError("Times must be in strictly increasing order.")

        self.rotations = rotations[:-1]
        self.rotvecs = (self.rotations.inv() * rotations[1:]).as_rotvec()

    def __call__(self, times: ArrayLike) -> Rotation:
        """Interpolate rotations.

        Compute the interpolated rotations at the given `times`.

        Parameters
        ----------
        times : array_like
            Times to compute the interpolations at. Can be a scalar or
            1-dimensional.

        Returns
        -------
        interpolated_rotation : `Rotation` instance
            Object containing the rotations computed at given `times`.

        """
        # Clearly differentiate from self.times property
        xp = array_namespace(self.times)
        compute_times = xp.asarray(
            times, device=device(self.times), dtype=self.times.dtype
        )
        if compute_times.ndim > 1:
            raise ValueError("`times` must be at most 1-dimensional.")

        single_time = compute_times.ndim == 0
        compute_times = xpx.atleast_nd(compute_times, ndim=1, xp=xp)

        # side = 'left' (default) excludes t_min.
        ind = xp.searchsorted(self.times, compute_times) - 1
        # Include t_min. Without this step, index for t_min equals -1
        ind = xp.where(compute_times == self.times[0], xp.asarray(0), ind)
        invalid_ind = xp.logical_or(ind < 0, ind > len(self.rotations) - 1)
        # We cannot error out on invalid indices for jit compiled code. To not produce an index
        # error, we set the index to 0 in case it is out of bounds, and later set the result to nan.
        # DECISION: Do we want to do this for all implementations or only for jax?
        if is_jax(xp):
            ind = xp.where(invalid_ind, xp.asarray(0), ind)
        else:
            if xp.any(invalid_ind):
                raise ValueError(
                    "Interpolation times must be within the range "
                    f"[{self.times[0]}, {self.times[-1]}], both inclusive."
                )
        alpha = (compute_times - self.times[ind]) / self.timedelta[ind]
        alpha = xp.where(invalid_ind, xp.asarray(xp.nan), alpha)

        result = self.rotations[ind] * Rotation.from_rotvec(
            self.rotvecs[ind[:, None], xp.arange(3)] * alpha[:, None]
        )

        if single_time:
            result = result[0]

        return result


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
        # Someone could have registered a different backend for jax, so we attempt to fetch the
        # updated backend here. If not, we fall back to the Array API backend.
        r._backend = backend_registry.get(array_namespace(c[0]), array_api_backend)
        r._quat = c[0]
        # We set _single to False for jax because the Array API backend supports broadcasting by
        # default and hence returns the correct shape without the _single workaround
        r._single = False
        return r

    register_pytree_node(Rotation, lambda v: ((v._quat,), None), rot_unflatten)
except ImportError:
    pass
