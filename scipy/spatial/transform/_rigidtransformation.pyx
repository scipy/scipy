# cython: cpow=True

import numpy as np
from scipy._lib._util import check_random_state, _transition_to_rng
from ._rotation import Rotation

cimport numpy as np
cimport cython

np.import_array()

cdef class RigidTransformation:
    """Rigid transformation in 3 dimensions.

    This class provides an interface to initialize from and represent rigid
    transformations (rotation and translation) in 3D space.

    The following operations on rigid transformations are supported:

    - Application on vectors
    - Transformation Composition
    - Transformation Inversion
    - Transformation Indexing

    Note that reflections are not supported, and coordinate systems must be
    right-handed. Because of this, this class more precisely represents
    *proper* rigid transformations [1]_ rather than rigid transformations more
    generally.

    Indexing within a transformation is supported since multiple transformations
    can be stored within a single `RigidTransformation` instance.

    To create `RigidTransformation` objects use ``from_...`` methods
    (see examples below). ``RigidTransformation(...)`` is not supposed
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
    from_rottrans
    from_expcoords
    from_dualquat
    as_matrix
    as_rottrans
    as_expcoords
    as_dualquat
    concatenate
    apply
    inv
    identity
    random

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Rigid_transformation

    Notes
    -----
    .. versionadded:: 1.16.0

    Examples
    --------

    >>> from scipy.spatial.transform import RigidTransformation as T
    >>> from scipy.spatial.transform import Rotation as R
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

    A `RigidTransformation` instance can be initialized in any of the
    above formats and converted to any of the others. The underlying object is
    independent of the representation used for initialization.

    A rigid transformation embodies the location and orientation of a coordinate
    frame in reference to another coordinate frame in 3D space.

    When we name transformations, we read the subscripts from right to left.
    So ``tAB`` represents a transformation A<-B and can be interpreted as:

    - the coordinates of B relative to A
    - the transformation of points from B to A
    - the pose of B described in A's coordinate system

        tAB
         ^^
         ||
         |--- from B
         |
         ---- to A

    When composing transformations, the order is important. Transformations
    are not commutative, so in general ``tAB * tBC`` is not the same as
    ``tBC * tAB``. Transformations are composed and applied to vectors
    right-to-left. So ``(tAB * tBC).apply(pC)`` is the same as
    ``tAB.apply(tBC.apply(pC))``.

    When composed, transformations should be ordered such that the
    multiplication operator is surrounded by a single frame, so that the frame
    "cancels out" and the outside frames are left. In the example below, B
    cancels out and the outside frames A and C are left. Or to put it another
    way, A <- C is the same as A <- B <- C.

                ----------- B cancels out
                |    |
                v    v
        tAC = tAB * tBC
               ^      ^
               |      |
               ------------ to A, from C are left

    The point pB defined in B is transformed to pA defined in A by:

               ------------ B cancels out
               |        |
               v        v
        pA = tAB.apply(pB)
              ^
              |
              ------------- A is left

    The following function can be used to plot transformations with Matplotlib
    by showing how they transform the standard x, y, z coordinate axes:

    >>> colors = ("#FF6666", "#005533", "#1199EE")  # Colorblind-safe RGB
    >>> def plot_transformed_axes(ax, t, name=None, scale=1):
    ...     r = t.rotation
    ...     d = t.translation
    ...     loc = np.array([d, d])
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
    ...         text_plot = text_loc_rot + d
    ...         ax.text(*text_plot, axlabel.upper(), color=c,
    ...                 va="center", ha="center")
    ...     ax.text(*t.translation, name, color="k", va="center", ha="center",
    ...             bbox={"fc": "w", "alpha": 0.8, "boxstyle": "circle"})

    Let's work through an example.

    First, define the "world frame" A, also called the "base frame".
    All frames are the identity transformation from their own perspective.

    >>> tA = T.identity()

    We will visualize a new frame B in A's coordinate system. So we need to
    define the transformation that converts coordinates from frame B to frame
    A (A<-B).

    Physically, let's imagine constructing B from A by:
    1) Rotating A by +90 degrees around its x-axis.
    2) Translating the rotated frame +2 units in A's -x direction.

    From A's perspective, B is at [-2, 0, 0] and rotated +90° about the x-axis,
    which is exactly the transform A<-B.

    >>> rAB = R.from_euler('xyz', [90, 0, 0], degrees=True)
    >>> dAB = np.array([-2, 0, 0])
    >>> tAB = T.from_rottrans(rAB, dAB)

    "tBA" is the inverse, i.e. the transformation B<-A.

    >>> tBA = tAB.inv()

    Let's plot these frames.

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tA, name="tA")    # A plotted in A
    >>> plot_transformed_axes(ax, tAB, name="tAB")  # B plotted in A
    >>> ax.set_title("A, B frames with respect to A")
    >>> ax.set_aspect("equal")
    >>> plt.show()

    Now let's visualize a new frame C in B's coordinate system.
    Let's imagine constructing C from B by:

    1) Translating B by 2 units in its +zB direction.
    2) Rotating B by 30 degrees around its z-axis.

    >>> dBC = np.array([0, 0, 2])
    >>> rBC = R.from_euler('xyz', [0, 0, 30], degrees=True)
    >>> tBC = T.from_rottrans(rBC, dBC)
    >>> tCB = tBC.inv()

    In order to plot these frames from a consistent perspective, we need to
    calculate the transformation between A and C. Note that we do not make this
    transformation directly, but instead compose intermediate transformations
    that let us get from A to C:

    >>> tCA = tCB * tBA  # C <- B <- A
    >>> tAC = tCA.inv()

    Now we can plot these three frames from A's perspective.

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tA, name="tA")    # A plotted in A
    >>> plot_transformed_axes(ax, tAB, name="tAB")  # B plotted in A
    >>> plot_transformed_axes(ax, tAC, name="tAC")  # C plotted in A
    >>> ax.set_title("A, B, C frames with respect to A")
    >>> ax.set_aspect("equal")
    >>> plt.show()

    Let's transform a point from A to B and C. We define a point in A and
    use the transformations we defined to transform it to coordinates in B and
    C:

    >>> p1_A = np.array([1, 0, 0])  # +1 in xA direction
    >>> p1_B = tBA.apply(p1_A)
    >>> p1_C = tCA.apply(p1_A)
    >>> print(p1_A)  # Original point in A
    [1 0 0]
    >>> print(p1_B)  # Point in B
    [3 0 0]
    >>> print(p1_C)  # Point in C
    [ 2.59807621 -1.5       -2.        ]

    We can also do the reverse. We define a point in C and transform it to A:

    >>> p2_C = np.array([0, 1, 0])  # +1 in yC direction
    >>> p2_A = tAC.apply(p2_C)
    >>> print(p2_C)  # Original point in C
    [0 1 0]
    >>> print(p2_A)  # Point in A
    [-2.5       -2.         0.8660254]

    Plot the frames with respect to A again, but also plot these two points:

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tA, name="tA")    # A plotted in A
    >>> plot_transformed_axes(ax, tAB, name="tAB")  # B plotted in A
    >>> plot_transformed_axes(ax, tAC, name="tAC")  # C plotted in A
    >>> ax.scatter(p1_A[0], p1_A[1], p1_A[2], color=colors[0])  # +1 xA
    >>> ax.scatter(p2_A[0], p2_A[1], p2_A[2], color=colors[1])  # +1 yC
    >>> ax.set_title("A, B, C frames and points with respect to A")
    >>> ax.set_aspect("equal")
    >>> plt.show()

    Up to this point, we have been visualizing frames from A's perspective.
    Let's use the transformations we defined to visualize the frames from C's
    perspective.

    Now C is the "base frame" or "world frame". All frames are the identity
    transformation from their own perspective.

    >>> tC = T.identity()

    We've already defined the transformations tCB and tCA, so we can plot
    everything from C's perspective:

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tC, name="tC")    # C plotted in C
    >>> plot_transformed_axes(ax, tCB, name="tCB")  # B plotted in C
    >>> plot_transformed_axes(ax, tCA, name="tCA")  # A plotted in C
    >>> ax.scatter(p1_C[0], p1_C[1], p1_C[2], color=colors[0])
    >>> ax.scatter(p2_C[0], p2_C[1], p2_C[2], color=colors[1])
    >>> ax.set_title("A, B, C frames and points with respect to C")
    >>> ax.set_aspect("equal")
    >>> plt.show()
    """

    cdef double[:, :, :] _matrix
    cdef bint _single

    def __init__(self, matrix, normalize=True, copy=True):
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

        if normalize or copy:
            matrix = matrix.copy()

        # Rigid transformations have the following matrix representation:
        # [R | t]
        # [0 | 1]
        # where R is a 3x3 orthonormal rotation matrix and t is a 3x1 translation
        # vector. The last row is always [0, 0, 0, 1] exactly

        # Check the last row
        last_row_ok = np.all(np.isclose(matrix[:, 3, :], [0, 0, 0, 1]), axis=1)

        # Check the determinant of the rotation matrix
        # (should be positive for right-handed rotations)
        dets = np.linalg.det(matrix[:, :3, :3])

        # Gramian orthonormality check
        # (should be the identity matrix for orthonormal matrices)
        grams = matrix[:, :3, :3] @ np.swapaxes(matrix[:, :3, :3], 1, 2)
        orthonormal_ok = np.all(np.isclose(grams, np.eye(3)), axis=(1, 2))

        for ind in range(num_transformations):
            if not last_row_ok[ind]:
                raise ValueError(f"Expected last row of transformation matrix {ind} to be "
                                 f"[0, 0, 0, 1], got {matrix[ind, 3, :]}.")

            # Check that the rotation matrix is orthonormal
            if not np.isclose(dets[ind], 1) or not orthonormal_ok[ind]:
                if normalize:
                    if dets[ind] <= 0:
                        raise ValueError("Non-positive determinant in rotation component of "
                                         f"matrix {ind}: {matrix[ind]}.")
                    else:
                        # Orthonormalize the rotation matrix
                        U, _, Vt = np.linalg.svd(matrix[ind, :3, :3])
                        matrix[ind, :3, :3] = U @ Vt
                else:
                    raise ValueError("Expected rotation component of transformation "
                                     f"matrix {ind} be orthonormal: {matrix[ind]}.")

        self._matrix = matrix

    @cython.embedsignature(True)
    @classmethod
    def from_rotation(cls, rotation):
        """Initialize from a rotation, without a translation.

        When applying this transformation to a vector, the result is the same as
        if the rotation was applied to the vector.
        ``t.from_rotation(r).apply(vector) == r.apply(vector)``

        Parameters
        ----------
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.

        Returns
        -------
        transformation : `RigidTransformation` instance

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating a transformation from a single rotation:

        >>> r = R.from_euler("ZYX", [90, 30, 0], degrees=True)
        >>> r.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> t = T.from_rotation(r)
        >>> t.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> t.single
        True

        The upper 3x3 submatrix of the transformation matrix is the rotation
        matrix:

        >>> np.allclose(t.as_matrix()[:3, :3], r.as_matrix(), atol=1e-12)
        True

        Creating multiple transformations from a stack of rotations:

        >>> r = R.from_euler("ZYX", [[90, 30, 0], [45, 30, 60]], degrees=True)
        >>> r.apply([1, 0, 0])
        array([[0.        , 0.8660254 , -0.5       ],
               [0.61237244, 0.61237244, -0.5       ]])
        >>> t = T.from_rotation(r)
        >>> t.apply([1, 0, 0])
        array([[0.        , 0.8660254 , -0.5       ],
               [0.61237244, 0.61237244, -0.5       ]])
        >>> t.single
        False
        >>> len(t)
        2
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
        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_translation(cls, translation):
        """Initialize from a translation numpy array, without a rotation.

        When applying this transformation to a vector, the result is the same as
        if the translation and vector were added together.
        ``t.from_translation(d).apply(vector) == d + vector``

        Parameters
        ----------
        translation : array_like, shape (N, 3) or (3,)
            A single translation vector or a stack of translation vectors.

        Returns
        -------
        transformation : `RigidTransformation` instance

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> import numpy as np

        Creating a transformation from a single translation vector:

        >>> d = np.array([2, 3, 4])
        >>> d + np.array([1, 0, 0])
        array([3, 3, 4])
        >>> t = T.from_translation(d)
        >>> t.apply([1, 0, 0])
        array([3., 3., 4.])
        >>> t.single
        True

        The rightmost column of the transformation matrix is the translation
        vector:

        >>> np.allclose(t.as_matrix()[:, 3], d)
        True

        Creating multiple transformations from a stack of translation vectors:

        >>> d = np.array([[2, 3, 4], [1, 0, 0]])
        >>> d + np.array([1, 0, 0])
        array([[3., 3., 4.],
               [2., 0., 0.]])
        >>> t = T.from_translation(d)
        >>> t.apply([1, 0, 0])
        array([[3., 3., 4.],
               [2., 0., 0.]])
        >>> t.single
        False
        >>> len(t)
        2
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
        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_matrix(cls, matrix):
        """Initialize from a 4x4 transformation matrix.

        Parameters
        ----------
        matrix : array_like, shape (4, 4) or (N, 4, 4)
            A single transformation matrix or a stack of transformation matrices.

        Returns
        -------
        transformation : `RigidTransformation` instance

        Notes
        -----
        4x4 rigid transformation matrices are of the form:

            [R | t]
            [0 | 1]

        where ``R`` is a 3x3 rotation matrix and ``t`` is a 3x1 translation vector
        ``[tx, ty, tz]``. As rotation matrices must be proper orthogonal, the
        rotation component is orthonormalized using singular value decomposition
        before initialization.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> import numpy as np

        Creating a transformation from a single matrix:

        >>> m = np.array([[0, 1, 0, 2],
        ...               [0, 0, 1, 3],
        ...               [1, 0, 0, 4],
        ...               [0, 0, 0, 1]])
        >>> t = T.from_matrix(m)
        >>> t.as_matrix()
        array([[0., 1., 0., 2.],
               [0., 0., 1., 3.],
               [1., 0., 0., 4.],
               [0., 0., 0., 1.]])
        >>> t.single
        True

        Creating a transformation from a stack of matrices:

        >>> m = np.array([np.eye(4), np.eye(4)])
        >>> t = T.from_matrix(m)
        >>> t.as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]],
               [[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        >>> t.single
        False
        >>> len(t)
        2

        Matrices with a rotation component that is not proper orthogonal are
        orthonormalized using singular value decomposition before
        initialization:

        >>> t = T.from_matrix(np.diag([2, 2, 2, 1]))
        >>> t.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])
        """
        return cls(matrix, normalize=True, copy=True)

    @cython.embedsignature(True)
    @classmethod
    def from_rottrans(cls, rotation, translation):
        """Initialize from a rotation and translation.

        When composing a rotation and translation, the rotation is applied
        before the translation, such that
        ``t = T.from_rottrans(rotation, translation)``
        is equivalent to
        ``t = T.from_rotation(rotation) * T.from_translation(translation)``.

        When applying a transformation to a vector, the result is the same as
        if the transformation was applied to the vector in the following way:
        ``t.apply(vector) == translation + rotation.apply(vector)``

        Parameters
        ----------
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.
        translation : array_like, shape (N, 3) or (3,)
            A single translation vector or a stack of translation vectors.

        Returns
        -------
        `RigidTransformation`, either a single transformation or a stack of
        transformations.
            - If rotation is single and translation is shape (3,), then a single
              transformation is returned.
            - Otherwise, a stack of transformations is returned.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating from a single rotation and translation:

        >>> r = R.from_euler("ZYX", [90, 30, 0], degrees=True)
        >>> r.as_matrix()
        array([[ 0.       , -1.,  0.        ],
               [ 0.8660254,  0.,  0.5       ],
               [-0.5      ,  0.,  0.8660254 ]])
        >>> d = np.array([2, 3, 4])
        >>> t = T.from_rottrans(r, d)
        >>> t.rotation.as_matrix()
        array([[ 0.       , -1.,  0.        ],
               [ 0.8660254,  0.,  0.5       ],
               [-0.5      ,  0.,  0.8660254 ]])
        >>> t.translation
        array([2., 3., 4.])
        >>> t.single
        True

        When applying a transformation to a vector, the result is the same as
        if the transformation was applied to the vector in the following way:
        ``t.apply(vector) == translation + rotation.apply(vector)``

        >>> r.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> d + r.apply([1, 0, 0])
        array([2.       , 3.8660254,  3.5     ])
        >>> t.apply([1, 0, 0])
        array([2.       , 3.8660254,  3.5     ])
        """
        return cls.from_translation(translation) * cls.from_rotation(rotation)

    @cython.embedsignature(True)
    @classmethod
    def from_expcoords(cls, expcoords):
        """Initialize from exponential coordinates.

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

        return cls(matrix, normalize=False, copy=False)

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
        transformation : `RigidTransformation` instance
            The identity transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> import numpy as np

        Creating a single identity transformation:

        >>> t = T.identity()
        >>> t.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])
        >>> t.single
        True

        The identity transformation can be applied to a vector without effect:

        >>> t.apply([1, 2, 3])
        array([1., 2., 3.])

        The identity transformation when composed with another transformation
        has no effect:

        >>> t = T.random()
        >>> np.allclose((T.identity() * t).as_matrix(),
        ...             t.as_matrix(), atol=1e-12)
        True

        Multiple identity transformations can be generated at once:

        >>> T.identity(2).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]],
               [[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        >>> t.single
        False
        >>> len(t)
        2
        """
        if num is None:  # single
            return cls(np.eye(4, dtype=float), normalize=False, copy=False)
        else:
            return cls(np.eye(4, dtype=float)[None, :, :].repeat(num, axis=0),
                       normalize=False, copy=False)

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
        transformation : `RigidTransformation` instance
            A single transformation or a stack of transformations.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> import numpy as np
        >>> rng = np.random.default_rng(seed=123)
        >>> t = T.random(rng=rng)
        >>> t.as_matrix()
        array([[-0.27687894,  0.08111092, -0.95747536,  0.21419513],
               [ 0.43673235, -0.87694466, -0.20058143,  0.19815685],
               [-0.85592225, -0.47369724,  0.20738374, -0.95648017],
               [ 0.        ,  0.        ,  0.        ,  1.        ]])
        >>> t.single
        True

        The translation component will be a random unit vector:

        >>> np.linalg.norm(t.translation)
        1.0

        Multiple transformations can be generated at once:

        >>> t = T.random(3)
        >>> t.single
        False
        >>> len(t)
        3
        """
        rng = check_random_state(rng)
        r = Rotation.random(num, rng=rng)
        t = Rotation.random(num, rng=rng).apply([1, 0, 0])
        return cls.from_rottrans(r, t)

    @cython.embedsignature(True)
    @classmethod
    def concatenate(cls, transformations):
        """
        Concatenate a sequence of `RigidTransformation` objects into a
        single object.

        Parameters
        ----------
        transformations : sequence of `RigidTransformation` objects, or
            a single `RigidTransformation` object.

        Returns
        -------
        transformation : `RigidTransformation` instance
            The concatenated transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> t1 = T.from_translation([1, 0, 0])
        >>> t2 = T.from_translation([[2, 0, 0], [3, 0, 0]])
        >>> T.concatenate([t1, t2]).translation
        array([[1., 0., 0.],
               [2., 0., 0.],
               [3., 0., 0.]])
        """
        if isinstance(transformations, RigidTransformation):
            return cls(transformations._matrix, normalize=False, copy=True)

        if not all(isinstance(x, RigidTransformation) for x in transformations):
            raise TypeError("input must contain RigidTransformation objects only")

        ms = [x.as_matrix()[np.newaxis, :, :] if x.single else x.as_matrix()
              for x in transformations]
        return cls(np.concatenate(ms), normalize=False, copy=False)

    @cython.embedsignature(True)
    def as_matrix(self):
        """Return a copy of the matrix representation of the transformation.

        4x4 rigid transformation matrices are of the form:

            [R | t]
            [0 | 1]

        where ``R`` is a 3x3 orthonormal rotation matrix and ``t`` is a 3x1
        translation vector ``[tx, ty, tz]``.

        Returns
        -------
        matrix : numpy.ndarray, shape (4, 4) or (N, 4, 4)
            A single transformation matrix or a stack of transformation matrices.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        A transformation matrix is a 4x4 matrix formed from a 3x3 rotation
        matrix and a 3x1 translation vector:

        >>> r = R.from_matrix([[0, 0, 1],
        ...                    [1, 0, 0],
        ...                    [0, 1, 0]])
        >>> d = np.array([2, 3, 4])
        >>> t = T.from_rottrans(r, d)
        >>> t.as_matrix()
        array([[ 0., 0., 1., 2.],
               [ 1., 0., 0., 3.],
               [ 0., 1., 0., 4.],
               [ 0., 0., 0., 1.]])

        >>> T.identity(2).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]],
               [[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        """
        if self._single:
            return np.array(self._matrix[0])
        else:
            return np.array(self._matrix)

    @cython.embedsignature(True)
    def as_rottrans(self):
        """Return the rotation and translation of the transformation, where
        the rotation is applied first, followed by the translation.

        4x4 rigid transformation matrices are of the form:

            [R | t]
            [0 | 1]

        Where ``R`` is a 3x3 orthonormal rotation matrix and ``t`` is a 3x1
        translation vector ``[tx, ty, tz]``. This function returns the rotation
        corresponding to the this rotation matrix ``Rotation.from_matrix(R)``
        and the translation vector ``t``.

        Returns
        -------
        rotation : `Rotation` instance
            The rotation of the transformation.
        translation : numpy.ndarray, shape (N, 3) or (3,)
            The translation of the transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Recover the rotation and translation from a transformation:

        >>> r = R.from_matrix([[0, 0, 1],
        ...                    [1, 0, 0],
        ...                    [0, 1, 0]])
        >>> d = np.array([2, 3, 4])
        >>> t = T.from_rottrans(r, d)
        >>> t_r, t_d = t.as_rottrans()
        >>> t_r.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])
        >>> t_d
        array([2., 3., 4.])

        The transformation applied to a vector is equivalent to the rotation
        applied to the vector followed by the translation:

        >>> r.apply([1, 0, 0])
        array([0., 1., 0.])
        >>> d + r.apply([1, 0, 0])
        array([2., 4., 4.])
        >>> t.apply([1, 0, 0])
        array([2., 4., 4.])
        """
        return self.rotation, self.translation

    @cython.embedsignature(True)
    def as_expcoords(self):
        """Return the exponential coordinates of the transformation.

        Returns
        -------
        expcoords : numpy.ndarray, shape (N, 6) or (6,)
            A single exponential coordinate vector or a stack of exponential
            coordinate vectors. The first three components define the
            rotation and the last three components define the translation.
        """
        expcoords = np.empty((self._matrix.shape[0], 6), dtype=float)
        rotations = Rotation.from_matrix(self._matrix[:, :3, :3])
        expcoords[:, :3] = rotations.as_rotvec()
        thetas = np.linalg.norm(expcoords[:, :3], axis=-1)
        nonzero_angle = thetas != 0.0
        expcoords[nonzero_angle, :3] /= thetas[nonzero_angle, np.newaxis]

        thetas = np.maximum(thetas, np.finfo(float).tiny)
        ti = 1.0 / thetas
        tan_term = -0.5 / np.tan(thetas / 2.0) + ti
        o0 = expcoords[:, 0]
        o1 = expcoords[:, 1]
        o2 = expcoords[:, 2]
        p0 = self._matrix[:, 0, 3]
        p1 = self._matrix[:, 1, 3]
        p2 = self._matrix[:, 2, 3]
        o00 = o0 * o0
        o01 = o0 * o1
        o02 = o0 * o2
        o11 = o1 * o1
        o12 = o1 * o2
        o22 = o2 * o2
        expcoords[:, 3] = (p0 * ((-o11 - o22) * tan_term + ti)
                           + p1 * (o01 * tan_term + 0.5 * o2)
                           + p2 * (o02 * tan_term - 0.5 * o1)
                           )
        expcoords[:, 4] = (p0 * (o01 * tan_term - 0.5 * o2)
                           + p1 * ((-o00 - o22) * tan_term + ti)
                           + p2 * (0.5 * o0 + o12 * tan_term)
                           )
        expcoords[:, 5] = (p0 * (o02 * tan_term + 0.5 * o1)
                           + p1 * (-0.5 * o0 + o12 * tan_term)
                           + p2 * ((-o00 - o11) * tan_term + ti)
                           )

        expcoords *= thetas[:, np.newaxis]

        # pure translations
        traces = np.einsum("nii", self._matrix[:, :3, :3])
        ind_only_translation = np.where(traces >= 3.0 - np.finfo(float).eps)[0]
        for i in ind_only_translation:
            expcoords[i, :3] = 0.0
            expcoords[i, 3:] = self._matrix[i, :3, 3]

        if self._single:
            return expcoords[0]
        else:
            return expcoords

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

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> t = T.random(3)
        >>> len(t)
        3

        >>> t = T.from_translation([1, 0, 0])
        >>> len(t)  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        TypeError: Single transformation has no len().
        """
        if self._single:
            raise TypeError("Single transformation has no len().")

        return self._matrix.shape[0]

    @cython.embedsignature(True)
    def __getitem__(self, indexer):
        """Extract transformation(s) at given index(es) from this object.

        Creates a new `RigidTransformation` instance containing a subset
        of transformations stored in this object.

        Parameters
        ----------
        indexer : int or slice or array_like
            Specifies which transformation(s) to extract. A single indexer must
            be specified, i.e. as if indexing a 1 dimensional array or list.

        Returns
        -------
        transformation : `RigidTransformation` instance
            Contains
                - a single transformation, if `indexer` is a single index
                - a stack of transformation(s), if `indexer` is a slice, or an
                  index array.

        Raises
        ------
        TypeError
            If the transformation is a single transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> d = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]  # 3 translations
        >>> t = T.from_translation(d)

        A single index returns a single transformation:

        >>> t[0].as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])

        A slice returns a stack of transformations:

        >>> t[1:3].translation
        array([[1., 0., 0.],
               [2., 0., 0.]])

        An index array returns a stack of transformations:

        >>> t[[0, 2]].translation
        array([[0., 0., 0.],
               [2., 0., 0.]])
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

        value : `RigidTransformation` instance
            The transformation(s) to set.

        Raises
        ------
        TypeError
            If the transformation is a single transformation.
        """
        if self._single:
            raise TypeError("Single transformation is not subscriptable.")

        if not isinstance(value, self.__class__):
            raise TypeError("value must be a RigidTransformation object")

        arr = np.asarray(self._matrix)
        arr[indexer] = value.as_matrix()
        self._matrix = arr

    @cython.embedsignature(True)
    def __mul__(RigidTransformation self, RigidTransformation other):
        """Compose this transformation with the other.

        If `p` and `q` are two transformations, then the composition of 'q
        followed by p' is equivalent to `p * q`. In terms of transformation
        matrices, the composition can be expressed as
        ``p.as_matrix() @ q.as_matrix()``.

        In terms of translations and rotations, the composition when applied to
        a vector is equivalent to
        ``p.translation + p.rotation.apply(q.translation)
        + (p.rotation * q.rotation).apply(vector)``.

        Parameters
        ----------
        other : `RigidTransformation` instance
            Object containing the transformations to be composed with this one.

        Returns
        -------
        `RigidTransformation` instance
            The composed transformation. This function supports composition of
            multiple transformations at a time. The following cases are
            possible:

            - Either ``p`` or ``q`` contains a single transformation. In this
              case `composition` contains the result of composing each
              transformation in the other object with the single transformation.
            - Both ``p`` and ``q`` contain ``N`` transformations. In this case
              each transformation ``p[i]`` is composed with the corresponding
              transformation ``q[i]`` and `output` contains ``N`` transformations.
        """
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
    def __pow__(RigidTransformation self, int n):
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
        `RigidTransformation` instance
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

        A rigid transformation is a composition of a rotation and a translation,
        where the rotation is applied first, followed by the translation. So
        the inverse transformation is equivalent to the inverse translation
        followed by the inverse transformation.

        Returns
        -------
        `RigidTransformation` instance
            The inverse of this transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> import numpy as np

        A transformation composed with its inverse results in an identity
        transformation:

        >>> t = T.random()
        >>> np.allclose((t.inv() * t).as_matrix(), np.eye(4), atol=1e-12)
        True

        >>> (t.inv() * t).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])

        The inverse rigid transformation is the same as the inverse translation
        followed by the inverse rotation:

        >>> r, d = t.as_rottrans()
        >>> r_inv = r.inv()  # inverse rotation
        >>> d_inv = -d  # inverse translation
        >>> t_r_inv = T.from_rotation(r_inv)
        >>> t_d_inv = T.from_translation(d_inv)
        >>> np.allclose((t_r_inv * t_d_inv).as_matrix(), t.inv().as_matrix(),
        ...             atol=1e-12)
        True
        >>> (t_r_inv * t_d_inv * t).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        """
        r_inv = self.rotation.inv()
        t_inv = -r_inv.apply(self.translation)
        return self.__class__.from_rottrans(r_inv, t_inv)

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

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        The rotation component is extracted from the transformation:

        >>> r = R.random(3)
        >>> d = np.array([1, 0, 0])
        >>> t = T.from_rottrans(r, d)
        >>> np.allclose(t.rotation.as_matrix(), r.as_matrix())
        True
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

        Returns
        -------
        translation : numpy.ndarray, shape (N, 3) or (3,)
            A single translation vector or a stack of translation vectors.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as T
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        The translation component is extracted from the transformation:

        >>> r = R.random()
        >>> d = np.array([[1, 0, 0], [2, 0, 0], [3, 0, 0]])
        >>> t = T.from_rottrans(r, d)
        >>> np.allclose(t.translation, d)
        True
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
