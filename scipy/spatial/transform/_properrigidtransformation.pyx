# cython: cpow=True

import numpy as np
from scipy._lib._util import check_random_state, _transition_to_rng
from ._rotation import Rotation

cimport numpy as np
cimport cython

np.import_array()

cdef class ProperRigidTransformation:
    """Proper rigid transformation in 3 dimensions.

    This class provides an interface to initialize from and represent proper
    rigid transformations (rotation and translation) in 3D space.

    The following operations on proper rigid transformations are supported:

    - Application on vectors
    - Transformation Composition
    - Transformation Inversion
    - Transformation Indexing

    Indexing within a transformation is supported since multiple transformations
    can be stored within a single `ProperRigidTransformation` instance.

    To create `ProperRigidTransformation` objects use ``from_...`` methods
    (see examples below). ``ProperRigidTransformation(...)`` is not supposed
    to be instantiated directly.

    Attributes
    ----------
    single
    rotation
    translation

    Methods
    -------
    __len__
    __getitem__
    __mul__
    __pow__
    from_rotation
    from_translation
    from_matrix
    from_transrot
    from_expcoords
    from_dualquat
    as_matrix
    as_expcoords
    as_dualquat
    concatenate
    apply
    inv
    identity
    random

    Notes
    -----
    .. versionadded:: 1.16.0

    Examples
    --------

    A `ProperRigidTransformation` instance can be initialized in any of the
    above formats and converted to any of the others. The underlying object is
    independent of the representation used for initialization.

    TODO: flesh this out

    >>> from scipy.spatial.transform import ProperRigidTransformation as T
    >>> from scipy.spatial.transform import Rotation as R
    >>> import numpy as np

    The following function can be used to plot transformations with Matplotlib by
    showing how they transform the standard x, y, z coordinate axes:

    >>> import matplotlib.pyplot as plt

    >>> def plot_rotated_axes(ax, t, name=None, scale=1):
    ...     colors = ("#FF6666", "#005533", "#1199EE")  # Colorblind-safe RGB
    ...     r = t.rotation
    ...     loc = np.array([t.translation, t.translation])
    ...     for i, (axis, c) in enumerate(zip((ax.xaxis, ax.yaxis, ax.zaxis),
    ...                                       colors)):
    ...         axlabel = axis.axis_name
    ...         axis.set_label_text(axlabel)
    ...         axis.label.set_color(c)
    ...         axis.line.set_color(c)
    ...         axis.set_tick_params(colors=c)
    ...         line = np.zeros((2, 3))
    ...         line[1, i] = scale
    ...         line_rot = r.apply(line)
    ...         line_plot = line_rot + loc
    ...         ax.plot(line_plot[:, 0], line_plot[:, 1], line_plot[:, 2], c)
    ...         text_loc = line[1]*1.2
    ...         text_loc_rot = r.apply(text_loc)
    ...         text_plot = text_loc_rot + loc[0]
    ...         ax.text(*text_plot, axlabel.upper(), color=c,
    ...                 va="center", ha="center")
    ...     ax.text(*t.translation, name, color="k", va="center", ha="center",
    ...             bbox={"fc": "w", "alpha": 0.8, "boxstyle": "circle"})

    Create two transformations - the identity, and a rotation followed by a
    translation:

    >>> t0 = T.identity()
    >>> t1 = T.from_transrot([3, 0, 0], R.from_euler("ZYX", [90, 30, 0], degrees=True))

    Add both transformations to a single plot:

    >>> ax = plt.figure().add_subplot(projection="3d", proj_type="ortho")
    >>> plot_rotated_axes(ax, t0, name="t0")
    >>> plot_rotated_axes(ax, t1, name="t1")
    >>> _ = ax.annotate(
    ...     "t0: Identity Transformation\\n"
    ...     "t1: Translation followed by Rotation",
    ...     xy=(0.6, 0.7), xycoords="axes fraction", ha="left"
    ... )
    >>> ax.set(xlim=(-1.25, 4.25), ylim=(-1.25, 1.25), zlim=(-1.25, 1.25))
    >>> ax.set(xticks=range(-1, 5), yticks=[-1, 0, 1], zticks=[-1, 0, 1])
    >>> ax.set_aspect("equal", adjustable="box")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.tight_layout()

    Show the plot:

    >>> plt.show()
    """

    cdef double[:, :, :] _matrix
    cdef bint _single

    def __init__(self, matrix, copy=True):
        """Initialize from a 4x4 transformation matrix."""
        self._single = False
        matrix = np.asarray(matrix, dtype=float)
        if (matrix.ndim not in [2, 3]
                or matrix.shape[0] == 0
                or matrix.shape[-1] != 4
                or matrix.shape[-2] != 4):
            raise ValueError("Expected `matrix` to have shape (4, 4), or (N, 4, 4), "
                             f"got {matrix.shape}.")

        # If a single transformation matrix is given, convert it to a
        # 2D 1 x 4 x 4 matrix but set self._single to True so that we can
        # return appropriate objects in the `to_...` methods
        if matrix.shape == (4, 4):
            matrix = matrix[None, :, :]
            self._single = True

        cdef Py_ssize_t num_transformations = matrix.shape[0]

        if copy:
            matrix = matrix.copy()

        # Proper rigid transformations have the following matrix representation:
        # [R | t]
        # [0 | 1]
        # where R is a 3x3 rotation matrix and t is a 3x1 translation vector.
        # The last row is always [0, 0, 0, 1] exactly
        for ind in range(num_transformations):
            if not np.allclose(matrix[ind, 3, :], np.array([0, 0, 0, 1])):
                raise ValueError("Expected last row of transformation matrix to be "
                                 f"[0, 0, 0, 1], got {matrix[ind, 3, :]}.")
            if not np.isclose(np.linalg.det(matrix[ind, :3, :3]), 1):
                # TODO: do we want to allow for non-unit determinants,
                # like Rotation.from_matrix?
                raise ValueError("Expected rotation matrix to have determinant 1, "
                                 f"got {np.linalg.det(matrix[ind, :3, :3])}.")

        self._matrix = matrix

    @cython.embedsignature(True)
    @classmethod
    def from_rotation(cls, rotation):
        """Initialize from a rotation, without a translation.

        Parameters
        ----------
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.

        Returns
        -------
        transformation : `ProperRigidTransformation` instance
        """
        rotmat = rotation.as_matrix()
        if rotation.single:
            rotmat = rotmat[None, :, :]
        num_transformations = len(rotmat)
        matrix = np.zeros((num_transformations, 4, 4), dtype=float)
        matrix[:, :3, :3] = rotmat
        matrix[:, 3, 3] = 1
        if rotation.single:
            matrix = matrix[0]
        return cls(matrix, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_translation(cls, translation):
        """Initialize from a translation numpy array, without a rotation.

        Parameters
        ----------
        translation : array_like, shape (N, 3) or (3,)
            A single translation vector or a stack of translation vectors.

        Returns
        -------
        transformation : `ProperRigidTransformation` instance
        """
        translation = np.asarray(translation, dtype=float)

        if (translation.ndim not in [1, 2]
                or translation.shape[0] == 0
                or translation.shape[len(translation.shape) - 1] != 3):
            raise ValueError("Expected `translation` to have shape (3,), or (N, 3), "
                             f"got {translation.shape}.")

        # If a single translation vector is given, convert it to a 2D 1 x 3 matrix
        single = False
        if translation.shape == (3,):
            translation = translation[None, :]
            single = True

        num_translations = translation.shape[0]

        matrix = np.repeat(np.eye(4, dtype=float)[None, :, :], num_translations, axis=0)
        matrix[:, :3, 3] = translation

        if single:
            matrix = matrix[0]
        return cls(matrix, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_matrix(cls, matrix, copy=True):
        """Initialize from a 4x4 transformation matrix.

        Parameters
        ----------
        matrix : array_like, shape (4, 4) or (N, 4, 4)
            A single transformation matrix or a stack of transformation matrices.

        Returns
        -------
        transformation : `ProperRigidTransformation` instance

        Notes
        -----
        4x4 matrices are expected to be of the form::

            [R | t]
            [0 | 1]

        where ``R`` is a 3x3 rotation matrix and ``t`` is a 3x1 translation vector
        ``[tx, ty, tz]``.
        """
        return cls(matrix, copy=copy)

    @cython.embedsignature(True)
    @classmethod
    def from_transrot(cls, translation, rotation):
        """Initialize from a translation and rotation.

        When composing a translation and rotation, the translation is applied
        after the rotation, such that
        ``T0 = T.from_transrot(translation, rotation)``
        is equivalent to
        ``T0 = T.from_translation(translation) * T.from_rotation(rotation)``.

        When applying a transformation to a vector, the result is the same as
        if the transformation was applied to the vector in the following way:
        ``T0.apply(vector) == translation + rotation.apply(vector)``

        Parameters
        ----------
        translation : array_like, shape (N, 3) or (3,)
            A single translation vector or a stack of translation vectors.
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.
        """
        return cls.from_translation(translation) * cls.from_rotation(rotation)

    @cython.embedsignature(True)
    @classmethod
    def from_expcoords(cls, expcoords):
        """Initialize from exponential coordinates.

        TODO: implement

        Parameters
        ----------
        expcoords : array_like, shape (N, 6) or (6,)
            A single exponential coordinate vector or a stack of exponential
            coordinate vectors. The first three components define the
            rotation and the last three components define the translation.
        """
        expcoords = np.asarray(expcoords, dtype=float)

        if (expcoords.ndim not in [1, 2] or expcoords.shape[0] == 0
                or expcoords.shape[-1] != 6):
            raise ValueError(
                "Expected `expcoords` to have shape (6,), or (N, 6), "
                f"got {expcoords.shape}.")

        single = expcoords.ndim == 1
        expcoords = np.atleast_2d(expcoords)

        rotations = Rotation.from_rotvec(expcoords[:, :3])

        theta = np.linalg.norm(expcoords[:, :3], axis=-1)
        ind_only_translation = theta == 0.0

        translations = np.empty((len(expcoords), 3), dtype=float)

        if not np.all(ind_only_translation):
            theta[ind_only_translation] = 1.0
            screw_axes = expcoords / theta[:, np.newaxis]

            tms = theta - np.sin(theta)
            cm1 = np.cos(theta) - 1.0
            o0 = screw_axes[:, 0]
            o1 = screw_axes[:, 1]
            o2 = screw_axes[:, 2]
            v0 = screw_axes[:, 3]
            v1 = screw_axes[:, 4]
            v2 = screw_axes[:, 5]
            o01tms = o0 * o1 * tms
            o12tms = o1 * o2 * tms
            o02tms = o0 * o2 * tms
            o0cm1 = o0 * cm1
            o1cm1 = o1 * cm1
            o2cm1 = o2 * cm1
            o00tms = o0 * o0 * tms
            o11tms = o1 * o1 * tms
            o22tms = o2 * o2 * tms
            translations[..., 0] = (-v0 * (o11tms + o22tms - theta)
                + v1 * (o01tms + o2cm1)
                + v2 * (o02tms - o1cm1))
            translations[..., 1] = (v0 * (o01tms - o2cm1)
                - v1 * (o00tms + o22tms - theta)
                + v2 * (o0cm1 + o12tms))
            translations[..., 2] = (v0 * (o02tms + o1cm1)
                - v1 * (o0cm1 - o12tms)
                - v2 * (o00tms + o11tms - theta))

        translations[ind_only_translation] = expcoords[ind_only_translation, 3:]

        # TODO alternatively, we could
        # return cls.from_translation(translations) * cls.from_rotation(rotations)

        matrix = np.empty((len(expcoords), 4, 4), dtype=float)
        matrix[:, :3, :3] = rotations.as_matrix()
        matrix[:, :3, 3] = translations
        matrix[:, 3, :3] = 0.0
        matrix[:, 3, 3] = 1.0
        if single:
            matrix = matrix[0]

        return cls(matrix, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_dualquat(cls, dualquat):
        """Initialize from a dual quaternion.

        TODO: implement

        Parameters
        ----------
        dualquat : array_like, shape (N, 8) or (8,)
            A single dual quaternion or a stack of dual quaternions.
        """
        raise NotImplementedError("from_dualquat not implemented")

    @cython.embedsignature(True)
    @classmethod
    def identity(cls, num=None):
        """Initialize an identity transformation.

        Composition with the identity transformation has no effect, and applying
        the identity transformation to a vector has no effect.

        Parameters
        ----------
        num : int, optional
            Number of identity transformations to generate. If None (default),
            then a single transformation is generated.

        Returns
        -------
        transformation : `ProperRigidTransformation` instance
            The identity transformation.
        """
        if num is None:  # single
            return cls(np.eye(4, dtype=float))
        else:
            return cls(np.eye(4, dtype=float)[None, :, :].repeat(num, axis=0))

    @cython.embedsignature(True)
    @classmethod
    @_transition_to_rng('random_state', position_num=2)
    def random(cls, num=None, rng=None):
        """Initialize a random transformation.

        This results in a random rotation combined with a random translation
        vector, randomly sampled from the unit sphere.

        Parameters
        ----------
        num : int, optional
            Number of random transformations to generate. If None (default), then a
            single transformation is generated.
        rng : `numpy.random.Generator`, optional
            Pseudorandom number generator state. When `rng` is None, a new
            `numpy.random.Generator` is created using entropy from the
            operating system. Types other than `numpy.random.Generator` are
            passed to `numpy.random.default_rng` to instantiate a `Generator`.

        Returns
        -------
        transformation : `ProperRigidTransformation` instance
            A single transformation or a stack of transformations.
        """
        rng = check_random_state(rng)
        r = Rotation.random(num, rng=rng)
        t = Rotation.random(num, rng=rng).apply([1, 0, 0])
        return cls.from_transrot(t, r)

    @cython.embedsignature(True)
    @classmethod
    def concatenate(cls, transformations):
        """
        Concatenate a sequence of `ProperRigidTransformation` objects into a
        single object.

        Parameters
        ----------
        transformations : sequence of `ProperRigidTransformation` objects, or
            a single `ProperRigidTransformation` object.

        Returns
        -------
        transformation : `ProperRigidTransformation` instance
            The concatenated transformation.
        """
        if isinstance(transformations, ProperRigidTransformation):
            return cls(transformations._matrix, copy=True)

        if not all(isinstance(x, ProperRigidTransformation) for x in transformations):
            raise TypeError("input must contain ProperRigidTransformation objects only")

        ms = [x.as_matrix()[np.newaxis, :, :] if x.single else x.as_matrix()
              for x in transformations]
        return cls(np.concatenate(ms), copy=False)

    @cython.embedsignature(True)
    def as_matrix(self):
        """Return a copy of the matrix representation of the transformation.

        4x4 matrices are of the form::

            [R | t]
            [0 | 1]

        where ``R`` is a 3x3 rotation matrix and ``t`` is a 3x1 translation vector
        ``[tx, ty, tz]``.

        Returns
        -------
        matrix : numpy.ndarray, shape (4, 4) or (N, 4, 4)
            A single transformation matrix or a stack of transformation matrices.
        """
        if self._single:
            return np.array(self._matrix[0])
        else:
            return np.array(self._matrix)

    @cython.embedsignature(True)
    def as_expcoords(self):
        """Return the exponential coordinates of the transformation.

        TODO: implement

        Returns
        -------
        expcoords : numpy.ndarray, shape TODO
            A single exponential coordinate vector or a stack of exponential
            coordinate vectors.
        """
        raise NotImplementedError("as_expcoords not implemented")

    @cython.embedsignature(True)
    def as_dualquat(self):
        """Return the dual quaternion representation of the transformation.

        TODO: implement

        Returns
        -------
        dualquat : numpy.ndarray, shape (N, 8) or (8,)
            A single dual quaternion vector or a stack of dual quaternion vectors.
        """
        raise NotImplementedError("as_dualquat not implemented")

    @cython.embedsignature(True)
    def __len__(self):
        """Return the number of transformations in this object.

        Multiple transformations can be stored in a single instance.

        Returns
        -------
        length : int
            The number of transformations in this object.

        Raises
        ------
        TypeError
            If the transformation is a single transformation.
        """
        if self._single:
            raise TypeError("Single transformation has no len().")

        return self._matrix.shape[0]

    @cython.embedsignature(True)
    def __getitem__(self, indexer):
        """Extract transformation(s) at given index(es) from this object.

        Creates a new `ProperRigidTransformation` instance containing a subset
        of transformations stored in this object.

        Parameters
        ----------
        indexer : int or slice or array_like
            Specifies which transformation(s) to extract. A single indexer must
            be specified, i.e. as if indexing a 1 dimensional array or list.

        Returns
        -------
        transformation : `ProperRigidTransformation` instance
            Contains
                - a single transformation, if `indexer` is a single index
                - a stack of transformation(s), if `indexer` is a slice, or an
                  index array.

        Raises
        ------
        TypeError
            If the transformation is a single transformation.
        """
        if self._single:
            raise TypeError("Single transformation is not subscriptable.")

        # Convert memoryview to numpy array before indexing:
        arr = np.asarray(self._matrix)
        return self.__class__(arr[indexer], copy=False)

    @cython.embedsignature(True)
    def __setitem__(self, indexer, value):
        """Set transformation(s) at given index(es) in this object.

        Parameters
        ----------
        indexer : int or slice or array_like
            Specifies which transformation(s) to replace. A single indexer must
            be specified, i.e. as if indexing a 1 dimensional array or list.

        value : `ProperRigidTransformation` instance
            The transformation(s) to set.

        Raises
        ------
        TypeError
            If the transformation is a single transformation.
        """
        if self._single:
            raise TypeError("Single transformation is not subscriptable.")

        if not isinstance(value, self.__class__):
            raise TypeError("value must be a ProperRigidTransformation object")

        arr = np.asarray(self._matrix)
        arr[indexer] = value.as_matrix()
        self._matrix = arr

    @cython.embedsignature(True)
    def __mul__(ProperRigidTransformation self, ProperRigidTransformation other):
        """Compose this transformation with the other."""
        len_self = len(self._matrix)
        len_other = len(other._matrix)
        if not(len_self == 1 or len_other == 1 or len_self == len_other):
            raise ValueError("Expected equal number of transformations in both "
                             f"or a single transformation in either object, "
                             f"got {len_self} transformations in first and "
                             f"{len_other} transformations in second object.")

        result = np.matmul(self._matrix, other._matrix)
        if self._single and other._single:
            result = result[0]
        return self.__class__(result, copy=False)

    @cython.embedsignature(True)
    def __pow__(ProperRigidTransformation self, int n):
        """Compose this transformation with itself `n` times.

        If `n` is negative, the inverse of the transformation is composed with
        itself `n` times. In other words, ``p ** -abs(n) == p.inv() ** abs(n)``.

        Non-integer values of `n` are not currently supported.

        Parameters
        ----------
        n : int
            The number of times to compose the transformation with itself.

        Returns
        -------
        `ProperRigidTransformation` instance
            If the input Rotation ``p`` contains ``N`` multiple rotations, then
            the output will contain ``N`` rotations where the ``i`` th rotation
            is equal to ``p[i] ** n``

        Notes
        -----
        There are three notable cases: if ``n == 1`` then the original rotation
        is returned, if ``n == 0`` then the identity rotation is returned, and
        if ``n == -1`` then ``p.inv()`` is returned.
        """
        # Exact short-cuts
        if n == 0:
            return self.__class__.identity(None if self._single else len(self._matrix))
        elif n == -1:
            return self.inv()
        elif n == 1:
            if self._single:
                return self.__class__(self._matrix[0], copy=True)
            else:
                return self.__class__(self._matrix, copy=True)
        else:
            matrix = self.as_matrix()
            if self._single:
                matrix = matrix[None, :, :]
            for i in range(len(matrix)):
                matrix[i] = np.linalg.matrix_power(matrix[i], abs(n))
            if self._single:
                matrix = matrix[0]

            T = self.__class__(matrix, copy=False)
            if n < 0:
                return T.inv()
            else:
                return T

    @cython.embedsignature(True)
    def inv(self):
        """Invert this transformation.

        Composition of a transformation with its inverse results in an identity
        transformation.

        Returns
        -------
        `ProperRigidTransformation` instance
            The inverse of this transformation.
        """
        r_inv = self.rotation.inv()
        t_inv = -r_inv.apply(self.translation)
        return self.__class__.from_transrot(t_inv, r_inv)

    @cython.embedsignature(True)
    def apply(self, vector, inverse=False):
        """Apply the transformation to a vector.

        If the original frame transforms to the final frame by this transformation,
        then its application to a vector can be seen in two ways:

            - As a projection of vector components expressed in the final frame
              to the original frame.
            - As the physical transformation of a vector being glued to the
              original frame as it transforms. In this case the vector components
              are expressed in the original frame before and after the
              transformation.

        In terms of rotation matrices and translation vectors, this application
        is the same as ``self.translation + self.rotation.as_matrix() @ vectors``.

        Parameters
        ----------
        vector : array_like, shape (N, 3) or (3,)
            A single vector or a stack of vectors.
        inverse : bool, optional
            If True, the inverse of the transformation is applied to the vector.

        Returns
        -------
        transformed_vector : numpy.ndarray, shape (N, 3) or (3,)
            The transformed vector(s).
            Shape depends on the following cases:

                - If object contains a single transformation (as opposed to a
                  stack with a single transformation) and a single vector is
                  specified with shape ``(3,)``, then `transformed_vector` has
                  shape ``(3,)``.
                - In all other cases, `transformed_vector` has shape ``(N, 3)``,
                  where ``N`` is either the number of transformations or vectors.
        """
        vector = np.asarray(vector, dtype=float)
        if (vector.ndim not in [1, 2]
                or vector.shape[-1] != 3
                or vector.shape[0] == 0):
            raise ValueError("Expected vector to have shape (N, 3), or (3,), "
                             f"got {vector.shape}.")

        vector_single = False
        if vector.ndim == 1:
            vector = vector[None, :]
            vector_single = True

        vector = np.hstack([vector, np.ones((vector.shape[0], 1))])

        if inverse:
            m = self.inv().as_matrix()
            if self._single:
                m = m[None, :, :]
        else:
            m = self._matrix

        # This einsum performs matrix multiplication of each of the (4, 4)
        # matrices in `m` with the (4,) vectors in `vector`, with proper
        # broadcasting for different dimensions of `m` and `vector`.
        res = np.einsum('ijk,ik->ij', m, vector)[:, :3]
        if self._single and vector_single:
            return res[0]
        else:
            return res

    @cython.embedsignature(True)
    @property
    def rotation(self):
        """Return the rotation component of the transformation.

        A transformation is a composition of a rotation and a translation, such
        that when applied to a vector, the vector is first rotated and then
        translated. This property returns the rotation part of the transformation.

        Returns
        -------
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.
        """
        if self._single:
            return Rotation.from_matrix(self._matrix[0, :3, :3])
        else:
            return Rotation.from_matrix(self._matrix[:, :3, :3])

    @cython.embedsignature(True)
    @property
    def translation(self):
        """Return the translation component of the transformation.

        A transformation is a composition of a rotation and a translation, such
        that when applied to a vector, the vector is first rotated and then
        translated. This property returns the translation part of the
        transformation.
        """
        if self._single:
            return np.array(self._matrix[0, :3, 3])
        else:
            return np.array(self._matrix[:, :3, 3])

    @cython.embedsignature(True)
    @property
    def single(self):
        """Whether this instance represents a single transformation.

        Single transformations are not subscriptable, and do not have a length.

        Returns
        -------
        single : bool
            True if this instance represents a single transformation, False
            otherwise.
        """
        return self._single
