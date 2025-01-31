# cython: cpow=True

import numpy as np
from ._rotation import Rotation, compose_quat

cimport numpy as np
cimport cython

np.import_array()


cdef _create_skew_matrix(vec):
    """Create skew-symmetric (aka cross-product) matrix for stack of vectors."""
    result = np.zeros((len(vec), 3, 3))
    result[:, 0, 1] = -vec[:, 2]
    result[:, 0, 2] = vec[:, 1]
    result[:, 1, 0] = vec[:, 2]
    result[:, 1, 2] = -vec[:, 0]
    result[:, 2, 0] = -vec[:, 1]
    result[:, 2, 1] = vec[:, 0]
    return result


cdef _compute_se3_exp_translation_transform(rot_vec):
    """Compute transform matrix from se3 translation part to SE3 translation.
    
    The transform matrix depends on the rotation vector.
    """
    angle = np.linalg.norm(rot_vec, axis=1)
    mask = angle < 1e-3

    k1 = np.empty(len(rot_vec))
    k1[mask] = 0.5 - angle[mask]**2 / 24 + angle[mask]**4 / 720
    k1[~mask] = (1.0 - np.cos(angle[~mask])) / angle[~mask]**2

    k2 = np.empty(len(rot_vec))
    k2[mask] = 1/6 - angle[mask]**2 / 120 + angle[mask]**4 / 5040
    k2[~mask] = (angle[~mask] - np.sin(angle[~mask])) / angle[~mask]**3

    s = _create_skew_matrix(rot_vec)

    return np.identity(3) + k1[:, None, None] * s + k2[:, None, None] * s @ s


cdef _compute_se3_log_translation_transform(rot_vec):
    """Compute transform matrix from SE3 translation to se3 translation part.
    
    It is inverse of `_compute_se3_exp_translation_transform` in a closed analytical
    form.
    """
    angle = np.linalg.norm(rot_vec, axis=1)
    mask = angle < 1e-3

    k = np.empty(len(rot_vec))
    k[mask] = 1/12 + angle[mask]**2 / 720 + angle[mask]**4 / 30240
    k[~mask] = (1 - 0.5 * angle[~mask] / np.tan(0.5 * angle[~mask])) / angle[~mask]**2

    s = _create_skew_matrix(rot_vec)

    return np.identity(3) - 0.5 * s + k[:, None, None] * s @ s


cdef _quaternion_conjugate(quat):
    conjugate = np.copy(quat)
    conjugate[:, :3] = -conjugate[:, :3]
    return conjugate


cdef _compose_dual_quaternion(real_part1, dual_part1, real_part2, dual_part2):
    prod_real = compose_quat(real_part1, real_part2)
    prod_dual = (np.asarray(compose_quat(real_part1, dual_part2))
                 + np.asarray(compose_quat(dual_part1, real_part2)))
    return prod_real, prod_dual


cdef _normalize_dual_quaternion(real_part, dual_part):
    """Ensure that the norm is 1 and that real and dual part are orthogonal."""
    real_part = np.copy(real_part)
    dual_part = np.copy(dual_part)

    # compute the dual quaternion product of the input and its
    # component-wise quaternion conjugate
    real_conjugate = _quaternion_conjugate(real_part)
    dual_conjugate = _quaternion_conjugate(dual_part)
    prod_real, prod_dual = _compose_dual_quaternion(
        real_part, dual_part, real_conjugate, dual_conjugate)

    # special case: invalid real quaternion
    prod_real_norm = np.linalg.norm(prod_real, axis=1)
    invalid_real_mask = prod_real_norm == 0.0
    real_part[invalid_real_mask, :4] = [0., 0., 0., 1.]
    prod_real_norm[invalid_real_mask] = 1.0

    # compute normalization factor
    real_inv_sqrt = 1.0 / prod_real_norm
    dual_inv_sqrt = -0.5 * prod_dual * real_inv_sqrt[:, np.newaxis] ** 3

    # normalize dual quaternion
    real_part = real_inv_sqrt[:, np.newaxis] * real_part
    dual_part = real_inv_sqrt[:, np.newaxis] * dual_part + np.asarray(
        compose_quat(dual_inv_sqrt, real_part))

    return real_part, dual_part


def _create_transformation_matrix(translations, rotation_matrices, single):
    """Create a matrix from translations and rotation matrices.

    Parameters
    ----------
    translations : array_like, shape (N, 3) or (3,)
        A stack of translation vectors.
    rotation_matrices : array_like, shape (N, 3, 3) or (3, 3)
        A stack of rotation matrices.
    single : bool
        Whether the output should be a single matrix or a stack of matrices.

    Returns
    -------
    matrix : numpy.ndarray, shape (N, 4, 4)
        A stack of transformation matrices.
    """
    if translations.ndim == 1:
        translations = translations[np.newaxis, ...]
    if rotation_matrices.ndim == 2:
        rotation_matrices = rotation_matrices[np.newaxis, ...]
    if not len(rotation_matrices) == len(translations):
        raise ValueError("The number of rotation matrices and translations must "
                         "be the same.")

    matrix = np.empty((len(translations), 4, 4), dtype=float)
    matrix[:, :3, :3] = rotation_matrices
    matrix[:, :3, 3] = translations
    matrix[:, 3, :3] = 0.0
    matrix[:, 3, 3] = 1.0
    if single:
        return matrix[0]
    else:
        return matrix


cdef class RigidTransformation:
    """Rigid transformation in 3 dimensions.

    This class provides an interface to initialize from and represent rigid
    transformations (rotation and translation) in 3D space. In different
    fields, this type of transformation may be referred to as "*pose*"
    (especially in robotics), "*extrinsic parameters*", or the "*model matrix*"
    (especially in computer graphics), but the core concept is the same: a
    rotation and translation describing the orientation of one 3D coordinate
    frame relative to another. Mathematically, these transformations belong to
    the Special Euclidean group SE(3), which encodes rotation (SO(3)) plus
    translation.

    The following operations on rigid transformations are supported:

    - Application on vectors
    - Transformation composition
    - Transformation inversion
    - Transformation indexing

    Note that coordinate systems must be right-handed. Because of this, this
    class more precisely represents *proper* rigid transformations in SE(3)
    rather than rigid transformations in E(3) more generally [1]_.

    Indexing within a transformation is supported since multiple
    transformations can be stored within a single `RigidTransformation`
    instance.

    To create `RigidTransformation` objects use ``from_...`` methods
    (see examples below). ``RigidTransformation(...)`` is not supposed
    to be instantiated directly.

    For rigorous introductions to rigid transformations, see [2]_, [3]_, and
    [4]_.

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
    from_matrix
    from_rotation
    from_translation
    from_components
    from_exp_coords
    from_dual_quat
    as_matrix
    as_components
    as_exp_coords
    as_dual_quat
    concatenate
    apply
    inv
    identity

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Rigid_transformation
    .. [2] https://motion.cs.illinois.edu/RoboticSystems/CoordinateTransformations.html
    .. [3] https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    .. [4] Kevin M. Lynch and Frank C. Park, "Modern Robotics: Mechanics,
           Planning, and Control" Chapter 3.3, 2017, Cambridge University Press.
           https://hades.mech.northwestern.edu/images/2/25/MR-v2.pdf#page=107.31
    .. [5] Paul Furgale, "Representing Robot Pose: The good, the bad, and the
           ugly", June 9, 2014.
           https://rpg.ifi.uzh.ch/docs/teaching/2024/FurgaleTutorial.pdf

    Notes
    -----
    .. versionadded:: 1.16.0

    Examples
    --------
    A `RigidTransformation` instance can be initialized in any of the
    above formats and converted to any of the others. The underlying object is
    independent of the representation used for initialization.

    **Notation Conventions and Composition**

    The notation here largely follows the convention defined in [5]_.
    When we name transformations, we read the subscripts from right to left.
    So ``tf_A_B`` represents a transformation A <- B and can be interpreted as:

    - the coordinates and orientation of B relative to A
    - the transformation of points from B to A
    - the pose of B described in A's coordinate system

    .. parsed-literal::
        :class: highlight-none

        tf_A_B
           ^ ^
           | |
           | --- from B
           |
           ----- to A

    When composing transformations, the order is important. Transformations
    are not commutative, so in general ``tf_A_B * tf_B_C`` is not the same as
    ``tf_B_C * tf_A_B``. Transformations are composed and applied to vectors
    right-to-left. So ``(tf_A_B * tf_B_C).apply(p_C)`` is the same as
    ``tf_A_B.apply(tf_B_C.apply(p_C))``.

    When composed, transformations should be ordered such that the
    multiplication operator is surrounded by a single frame, so the frame
    "cancels out" and the outside frames A and C are left. In the example
    below, B cancels out and the outside frames A and C are left. Or to put
    it another way, A <- C is the same as A <- B <- C.

    .. parsed-literal::
        :class: highlight-none

                      ----------- B cancels out
                      |      |
                      v      v
        tf_A_C = tf_A_B * tf_B_C
                    ^          ^
                    |          |
                    ------------ to A, from C are left

    When we notate vectors, we write the subscript of the frame that the vector
    is defined in. So ``p_B`` means the point ``p`` defined in frame B. To
    transform this point from frame B to coordinates in frame A, we apply the
    transformation ``tf_A_B`` to the vector, lining things up such that the
    notated B frames are next to each other and "cancel out".

    .. parsed-literal::
        :class: highlight-none

                   ------------ B cancels out
                   |         |
                   v         v
        p_A = tf_A_B.apply(p_B)
                 ^
                 |
                 -------------- A is left

    **Visualization**

    >>> from scipy.spatial.transform import RigidTransformation as Tf
    >>> from scipy.spatial.transform import Rotation as R
    >>> import numpy as np

    The following function can be used to plot transformations with Matplotlib
    by showing how they transform the standard x, y, z coordinate axes:

    >>> import matplotlib.pyplot as plt
    >>> colors = ("#FF6666", "#005533", "#1199EE")  # Colorblind-safe RGB
    >>> def plot_transformed_axes(ax, tf, name=None, scale=1):
    ...     r = tf.rotation
    ...     t = tf.translation
    ...     loc = np.array([t, t])
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
    ...         text_plot = text_loc_rot + t
    ...         ax.text(*text_plot, axlabel.upper(), color=c,
    ...                 va="center", ha="center")
    ...     ax.text(*tf.translation, name, color="k", va="center", ha="center",
    ...             bbox={"fc": "w", "alpha": 0.8, "boxstyle": "circle"})

    **Defining Frames**

    Let's work through an example.

    First, define the "world frame" A, also called the "base frame".
    All frames are the identity transformation from their own perspective.

    >>> tf_A = Tf.identity()

    We will visualize a new frame B in A's coordinate system. So we need to
    define the transformation that converts coordinates from frame B to frame
    A (A <- B).

    Physically, let's imagine constructing B from A by:

    1) Rotating A by +90 degrees around its x-axis.
    2) Translating the rotated frame +2 units in A's -x direction.

    From A's perspective, B is at [-2, 0, 0] and rotated +90 degrees about the
    x-axis, which is exactly the transform A <- B.

    >>> t_A_B = np.array([-2, 0, 0])
    >>> r_A_B = R.from_euler('xyz', [90, 0, 0], degrees=True)
    >>> tf_A_B = Tf.from_components(t_A_B, r_A_B)

    Let's plot these frames.

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_A, name="tfA")     # A plotted in A
    >>> plot_transformed_axes(ax, tf_A_B, name="tfAB")  # B plotted in A
    >>> ax.set_title("A, B frames with respect to A")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()

    Now let's visualize a new frame C in B's coordinate system.
    Let's imagine constructing C from B by:

    1) Translating B by 2 units in its +z direction.
    2) Rotating B by +30 degrees around its z-axis.

    >>> t_B_C = np.array([0, 0, 2])
    >>> r_B_C = R.from_euler('xyz', [0, 0, 30], degrees=True)
    >>> tf_B_C = Tf.from_components(t_B_C, r_B_C)

    In order to plot these frames from a consistent perspective, we need to
    calculate the transformation between A and C. Note that we do not make this
    transformation directly, but instead compose intermediate transformations
    that let us get from C to A:

    >>> tf_A_C = tf_A_B * tf_B_C  # A <- B <- C

    Now we can plot these three frames from A's perspective.

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_A, name="tfA")     # A plotted in A
    >>> plot_transformed_axes(ax, tf_A_B, name="tfAB")  # B plotted in A
    >>> plot_transformed_axes(ax, tf_A_C, name="tfAC")  # C plotted in A
    >>> ax.set_title("A, B, C frames with respect to A")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()

    **Transforming Vectors**

    Let's transform a vector from A, to B and C. To do this, we will first
    invert the transformations we already have from B and C, to A.

    >>> tf_B_A = tf_A_B.inv()  # B <- A
    >>> tf_C_A = tf_A_C.inv()  # C <- A

    Now we can define a point in A and use the above transformations to get its
    coordinates in B and C:

    >>> p1_A = np.array([1, 0, 0])  # +1 in x_A direction
    >>> p1_B = tf_B_A.apply(p1_A)
    >>> p1_C = tf_C_A.apply(p1_A)
    >>> print(p1_A)  # Original point 1 in A
    [1 0 0]
    >>> print(p1_B)  # Point 1 in B
    [3. 0. 0.]
    >>> print(p1_C)  # Point 1 in C
    [ 2.59807621 -1.5       -2.        ]

    We can also do the reverse. We define a point in C and transform it to A:

    >>> p2_C = np.array([0, 1, 0])  # +1 in y_C direction
    >>> p2_A = tf_A_C.apply(p2_C)
    >>> print(p2_C)  # Original point 2 in C
    [0 1 0]
    >>> print(p2_A)  # Point 2 in A
    [-2.5       -2.         0.8660254]

    Plot the frames with respect to A again, but also plot these two points:

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_A, name="tfA")     # A plotted in A
    >>> plot_transformed_axes(ax, tf_A_B, name="tfAB")  # B plotted in A
    >>> plot_transformed_axes(ax, tf_A_C, name="tfAC")  # C plotted in A
    >>> ax.scatter(p1_A[0], p1_A[1], p1_A[2], color=colors[0])  # +1 x_A
    >>> ax.scatter(p2_A[0], p2_A[1], p2_A[2], color=colors[1])  # +1 y_C
    >>> ax.set_title("A, B, C frames and points with respect to A")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()

    **Switching Base Frames**

    Up to this point, we have been visualizing frames from A's perspective.
    Let's use the transformations we defined to visualize the frames from C's
    perspective.

    Now C is the "base frame" or "world frame". All frames are the identity
    transformation from their own perspective.

    >>> tf_C = Tf.identity()

    We've already defined the transformation C <- A, and can obtain C <- B by
    inverting the existing transformation B <- C.

    >>> tf_C_B = tf_B_C.inv()  # C <- B

    This lets us plot everything from C's perspective:

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_C, name="tfC")     # C plotted in C
    >>> plot_transformed_axes(ax, tf_C_B, name="tfCB")  # B plotted in C
    >>> plot_transformed_axes(ax, tf_C_A, name="tfCA")  # A plotted in C
    >>> ax.scatter(p1_C[0], p1_C[1], p1_C[2], color=colors[0])
    >>> ax.scatter(p2_C[0], p2_C[1], p2_C[2], color=colors[1])
    >>> ax.set_title("A, B, C frames and points with respect to C")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()
    """

    cdef double[:, :, :] _matrix
    cdef bint _single

    def __init__(self, matrix, normalize=True, copy=True):
        """Initialize from a 4x4 transformation matrix.

        Rotations are not meant to be initialized directly. Please use one of
        the `from_...` methods instead.

        Parameters
        ----------
        matrix : array_like, shape (4, 4) or (N, 4, 4)
            A single transformation matrix or a stack of transformation
            matrices.
        normalize : bool, optional
            If True, orthonormalize the rotation matrix using singular value
            decomposition.
        copy : bool, optional
            If True, copy the input matrix. If False, a reference to the input
            matrix is used. If normalize is True, the input matrix is always
            copied regardless of the value of copy.

        Returns
        -------
        transformation : `RigidTransformation` instance
        """
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

        if normalize or copy:
            matrix = matrix.copy()

        # Rigid transformations have the following matrix representation:
        # [R | t]
        # [0 | 1]
        # where R is a 3x3 orthonormal rotation matrix and t is a 3x1 translation
        # vector. The last row is always [0, 0, 0, 1] exactly

        # Check the last row. Exact match is required, as this last row should
        # not accumulate floating point error during composition.
        last_row_ok = np.all(matrix[:, 3, :] == np.array([0, 0, 0, 1]), axis=1)
        if np.any(~last_row_ok):
            ind = np.where(~last_row_ok)[0][0]
            raise ValueError(f"Expected last row of transformation matrix {ind} to be "
                             f"exactly [0, 0, 0, 1], got {matrix[ind, 3, :]}.")

        # The Rotation.from_matrix() method orthogonalizes the rotation
        # component of the transformation matrix. While this does have some
        # overhead in converting a rotation matrix to a quaternion and back, it
        # allows for skipping singular value decomposition for near-orthogonal
        # matrices, which is a computationally expensive operation.
        if normalize:
            matrix[:, :3, :3] = Rotation.from_matrix(matrix[:, :3, :3]).as_matrix()
        else:
            # Check for right-handedness
            dets = np.linalg.det(matrix[:, :3, :3])
            if np.any(dets < 0):
                idx = np.where(dets < 0)[0][0]
                raise ValueError("Non-positive determinant (left-handed or null "
                                 f"coordinate frame) in rotation component of "
                                 f"transformation matrix {idx}: "
                                 f"{matrix[idx, :3, :3]}.")

            # Gramian orthogonality check
            gramians = matrix[:, :3, :3] @ matrix[:, :3, :3].transpose(0, 2, 1)
            is_orthogonal = np.all(np.isclose(gramians, np.eye(3), atol=1e-12),
                                   axis=(1, 2))
            if np.any(~is_orthogonal):
                idx = np.where(~is_orthogonal)[0][0]
                raise ValueError("Rotation component of transformation matrix "
                                 f"{idx} is not orthogonal: {matrix[idx, :3, :3]}")

        self._matrix = matrix

    def __repr__(self):
        m = f"{self.as_matrix()!r}".splitlines()
        # bump indent (+32 characters)
        m[1:] = [" " * 32 + m[i] for i in range(1, len(m))]
        return "RigidTransformation.from_matrix(" + "\n".join(m) + ")"

    @cython.embedsignature(True)
    @classmethod
    def from_matrix(cls, matrix):
        """Initialize from a 4x4 transformation matrix.

        Parameters
        ----------
        matrix : array_like, shape (4, 4) or (N, 4, 4)
            A single transformation matrix or a stack of transformation
            matrices.

        Returns
        -------
        transformation : `RigidTransformation` instance

        Notes
        -----
        4x4 rigid transformation matrices are of the form:

        ..

            [R | t]
            [0 | 1]

        where ``R`` is a 3x3 rotation matrix and ``t`` is a 3x1 translation
        vector ``[tx, ty, tz]``. As rotation matrices must be proper
        orthogonal, the rotation component is orthonormalized using singular
        value decomposition before initialization.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> import numpy as np

        Creating a transformation from a single matrix:

        >>> m = np.array([[0, 1, 0, 2],
        ...               [0, 0, 1, 3],
        ...               [1, 0, 0, 4],
        ...               [0, 0, 0, 1]])
        >>> tf = Tf.from_matrix(m)
        >>> tf.as_matrix()
        array([[0., 1., 0., 2.],
               [0., 0., 1., 3.],
               [1., 0., 0., 4.],
               [0., 0., 0., 1.]])
        >>> tf.single
        True

        Creating a transformation from a stack of matrices:

        >>> m = np.array([np.eye(4), np.eye(4)])
        >>> tf = Tf.from_matrix(m)
        >>> tf.as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]],
               [[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        >>> tf.single
        False
        >>> len(tf)
        2

        Matrices with a rotation component that is not proper orthogonal are
        orthonormalized using singular value decomposition before
        initialization:

        >>> tf = Tf.from_matrix(np.diag([2, 2, 2, 1]))
        >>> tf.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])
        """
        return cls(matrix, normalize=True, copy=True)

    @cython.embedsignature(True)
    @classmethod
    def from_rotation(cls, rotation):
        """Initialize from a rotation, without a translation.

        When applying this transformation to a vector ``v``, the result is the
        same as if the rotation was applied to the vector.
        ``Tf.from_rotation(r).apply(v) == r.apply(v)``

        Parameters
        ----------
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.

        Returns
        -------
        transformation : `RigidTransformation` instance

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating a transformation from a single rotation:

        >>> r = R.from_euler("ZYX", [90, 30, 0], degrees=True)
        >>> r.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> tf = Tf.from_rotation(r)
        >>> tf.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> tf.single
        True

        The upper 3x3 submatrix of the transformation matrix is the rotation
        matrix:

        >>> np.allclose(tf.as_matrix()[:3, :3], r.as_matrix(), atol=1e-12)
        True

        Creating multiple transformations from a stack of rotations:

        >>> r = R.from_euler("ZYX", [[90, 30, 0], [45, 30, 60]], degrees=True)
        >>> r.apply([1, 0, 0])
        array([[0.        , 0.8660254 , -0.5       ],
               [0.61237244, 0.61237244, -0.5       ]])
        >>> tf = Tf.from_rotation(r)
        >>> tf.apply([1, 0, 0])
        array([[0.        , 0.8660254 , -0.5       ],
               [0.61237244, 0.61237244, -0.5       ]])
        >>> tf.single
        False
        >>> len(tf)
        2
        """
        if not isinstance(rotation, Rotation):
            raise ValueError("Expected `rotation` to be a `Rotation` instance, "
                             f"got {type(rotation)}.")
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

        When applying this transformation to a vector ``v``, the result is the
        same as if the translation and vector were added together. If ``t`` is
        the displacement vector of the translation, then:

        ``Tf.from_translation(t).apply(v) == t + v``

        Parameters
        ----------
        translation : array_like, shape (N, 3) or (3,)
            A single translation vector or a stack of translation vectors.

        Returns
        -------
        transformation : `RigidTransformation` instance

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> import numpy as np

        Creating a transformation from a single translation vector:

        >>> t = np.array([2, 3, 4])
        >>> t + np.array([1, 0, 0])
        array([3, 3, 4])
        >>> tf = Tf.from_translation(t)
        >>> tf.apply([1, 0, 0])
        array([3., 3., 4.])
        >>> tf.single
        True

        The top 3x1 points in the rightmost column of the transformation matrix
        is the translation vector:

        >>> tf.as_matrix()
        array([[1., 0., 0., 2.],
               [0., 1., 0., 3.],
               [0., 0., 1., 4.],
               [0., 0., 0., 1.]])
        >>> np.allclose(tf.as_matrix()[:3, 3], t)
        True

        Creating multiple transformations from a stack of translation vectors:

        >>> t = np.array([[2, 3, 4], [1, 0, 0]])
        >>> t + np.array([1, 0, 0])
        array([[3, 3, 4],
               [2, 0, 0]])
        >>> tf = Tf.from_translation(t)
        >>> tf.apply([1, 0, 0])
        array([[3., 3., 4.],
               [2., 0., 0.]])
        >>> np.allclose(tf.as_matrix()[:, :3, 3], t)
        True
        >>> tf.single
        False
        >>> len(tf)
        2
        """
        translation = np.asarray(translation, dtype=float)

        if (translation.ndim not in [1, 2]
                or translation.shape[0] == 0
                or translation.shape[-1] != 3):
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
    def from_components(cls, translation, rotation):
        """Initialize a rigid transformation from translation and rotation
        components.

        When creating a rigid transformation from a translation and rotation,
        the translation is applied after the rotation, such that
        ``tf = Tf.from_components(translation, rotation)``
        is equivalent to
        ``tf = Tf.from_translation(translation) * Tf.from_rotation(rotation)``.

        When applying a transformation to a vector ``v``, the result is the
        same as if the transformation was applied to the vector in the
        following way: ``tf.apply(v) == translation + rotation.apply(v)``

        Parameters
        ----------
        translation : array_like, shape (N, 3) or (3,)
            A single translation vector or a stack of translation vectors.
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.

        Returns
        -------
        `RigidTransformation`
            If rotation is single and translation is shape (3,), then a single
            transformation is returned.
            Otherwise, a stack of transformations is returned.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating from a single rotation and translation:

        >>> t = np.array([2, 3, 4])
        >>> r = R.from_euler("ZYX", [90, 30, 0], degrees=True)
        >>> r.as_matrix()
        array([[ 0.       , -1.,  0.        ],
               [ 0.8660254,  0.,  0.5       ],
               [-0.5      ,  0.,  0.8660254 ]])
        >>> tf = Tf.from_components(t, r)
        >>> tf.rotation.as_matrix()
        array([[ 0.       , -1.,  0.        ],
               [ 0.8660254,  0.,  0.5       ],
               [-0.5      ,  0.,  0.8660254 ]])
        >>> tf.translation
        array([2., 3., 4.])
        >>> tf.single
        True

        When applying a transformation to a vector ``v``, the result is the
        same as if the transformation was applied to the vector in the
        following way: ``tf.apply(v) == translation + rotation.apply(v)``

        >>> r.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> t + r.apply([1, 0, 0])
        array([2.       , 3.8660254,  3.5     ])
        >>> tf.apply([1, 0, 0])
        array([2.       , 3.8660254,  3.5     ])
        """
        return cls.from_translation(translation) * cls.from_rotation(rotation)

    @cython.embedsignature(True)
    @classmethod
    def from_exp_coords(cls, exp_coords):
        r"""Initialize from exponential coordinates of transformation.

        This implements the exponential map that converts 6-dimensional real
        vectors to SE(3). The first three components of the vector encode
        the rotation and the last three components encode the translation.

        The exponential coordinates of transformation can be split up into
        a screw axis :math:`\mathcal{S} \in \mathbb{R}^6` and a scalar
        :math:`\theta` that defines the magnitude of the displacement. For
        pure translations, the norm of the last three components defines
        :math:`\theta` and for all other transformations, the norm of the first
        three components defines :math:`\theta`.

        Parameters
        ----------
        exp_coords : array_like, shape (N, 6) or (6,)
            A single exponential coordinate vector or a stack of exponential
            coordinate vectors. The first three components define the
            rotation and the last three components define the translation.

        Returns
        -------
        transformation : `RigidTransformation` instance
            A single transformation or a stack of transformations.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> import numpy as np

        Creating from a single 6d vector of exponential coordinates:

        >>> tf = Tf.from_exp_coords([
        ...     -2.01041204, -0.52983629, 0.65773501,
        ...     0.10386614, 0.05855009, 0.54959179])
        >>> tf.as_matrix()
        array([[0.76406621, 0.10504613, -0.63652819, -0.10209961],
               [0.59956454, -0.47987325, 0.64050295, 0.40158789],
               [-0.2381705, -0.87102639, -0.42963687, 0.19637636],
               [0., 0., 0., 1.]])
        >>> tf.single
        True

        A vector of zeros represents no transformation:

        >>> tf = Tf.from_exp_coords(np.zeros(6))
        >>> tf.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])

        The last three numbers encode translation. If the first three
        numbers are zero, the last three components can be interpreted as the
        translation:

        >>> tf_trans = Tf.from_exp_coords([0, 0, 0, 4.3, -2, 3.4])
        >>> tf_trans.translation
        array([4.3, -2., 3.4])

        The first three numbers encode rotation as a rotation vector:

        >>> tf_rot = Tf.from_exp_coords([0.5, 0.3, 0.1, 0, 0, 0])
        >>> tf_rot.as_matrix()
        array([[ 0.95144143, -0.02143004,  0.307083  ,  0.        ],
               [ 0.16710577,  0.87374771, -0.45677194,  0.        ],
               [-0.25852442,  0.48590709,  0.83490085,  0.        ],
               [ 0.        ,  0.        ,  0.        ,  1.        ]])

        >>> tf_rot.rotation.as_rotvec()
        array([0.5, 0.3, 0.1])

        Combining translation and rotation preserves the rotation vector,
        but changes the last three components as they encode translation and
        rotation:

        >>> (tf_trans * tf_rot).as_exp_coords()
        array([0.5, 0.3, 0.1, 3.64305882, -1.25879559, 4.46109265])
        """
        exp_coords = np.asarray(exp_coords, dtype=float)

        if (exp_coords.ndim not in [1, 2] or exp_coords.shape[0] == 0
                or exp_coords.shape[-1] != 6):
            raise ValueError(
                "Expected `exp_coords` to have shape (6,), or (N, 6), "
                f"got {exp_coords.shape}.")

        single = exp_coords.ndim == 1
        exp_coords = np.atleast_2d(exp_coords)

        rot_vec = exp_coords[:, :3]
        rotations = Rotation.from_rotvec(rot_vec)

        translations = np.einsum('ijk,ik->ij',
                                 _compute_se3_exp_translation_transform(rot_vec),
                                 exp_coords[:, 3:])
        matrix = _create_transformation_matrix(translations, rotations.as_matrix(),
                                               single)

        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_dual_quat(cls, dual_quat, *, scalar_first=False):
        """Initialize from a unit dual quaternion.

        Unit dual quaternions encode orientation in a real unit quaternion
        and translation in a dual quaternion. There is a double cover, i.e.,
        the unit dual quaternions q and -q represent the same transformation.

        Unit dual quaternions must have a real quaternion with unit norm and
        a dual quaternion that is orthogonal to the real quaternion to satisfy
        the unit norm constraint. This function will enforce both properties
        through normalization.

        Parameters
        ----------
        dual_quat : array_like, shape (N, 8) or (8,)
            A single unit dual quaternion or a stack of unit dual quaternions.
            The real part is stored in the first four components and the dual
            part in the last four components.
        scalar_first : bool, optional
            Whether the scalar component goes first or last in the two
            individual quaternions that represent the real and the dual part.
            Default is False, i.e. the scalar-last order is used.

        Returns
        -------
        transformation : `RigidTransformation` instance
            A single transformation or a stack of transformations.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> import numpy as np

        Creating from a single unit dual quaternion:

        >>> tf = Tf.from_dual_quat([
        ...     0.0617101, -0.06483886, 0.31432811, 0.94508498,
        ...     0.04985168, -0.26119618, 0.1691491, -0.07743254])
        >>> tf.as_matrix()
        array([[0.79398752, -0.60213598, -0.08376202, 0.24605262],
               [0.58613113, 0.79477941, -0.15740392, -0.4932833],
               [0.16135089, 0.07588122, 0.98397557, 0.34262676],
               [0., 0., 0., 1.]])
        >>> tf.single
        True
        """
        dual_quat = np.asarray(dual_quat, dtype=float)

        if (dual_quat.ndim not in [1, 2] or dual_quat.shape[0] == 0
                or dual_quat.shape[-1] != 8):
            raise ValueError(
                "Expected `dual_quat` to have shape (8,), or (N, 8), "
                f"got {dual_quat.shape}.")

        single = dual_quat.ndim == 1
        dual_quat = np.atleast_2d(dual_quat)

        real_part = dual_quat[:, :4]
        dual_part = dual_quat[:, 4:]
        if scalar_first:
            real_part = np.roll(real_part, -1, axis=1)
            dual_part = np.roll(dual_part, -1, axis=1)

        real_part, dual_part = _normalize_dual_quaternion(real_part, dual_part)

        rotation = Rotation.from_quat(real_part)
        translation = 2.0 * np.asarray(
            compose_quat(dual_part, rotation.inv().as_quat()))[:, :3]
        matrix = _create_transformation_matrix(translation, rotation.as_matrix(),
                                               single)

        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def identity(cls, num=None):
        """Initialize an identity transformation.

        Composition with the identity transformation has no effect, and
        applying the identity transformation to a vector has no effect.

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
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating a single identity transformation:

        >>> tf = Tf.identity()
        >>> tf.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])
        >>> tf.single
        True

        The identity transformation can be applied to a vector without effect:

        >>> tf.apply([1, 2, 3])
        array([1., 2., 3.])

        The identity transformation when composed with another transformation
        has no effect:

        >>> rng = np.random.default_rng(123)
        >>> t = rng.random(3)
        >>> r = R.random(rng=rng)
        >>> tf = Tf.from_components(t, r)
        >>> np.allclose((Tf.identity() * tf).as_matrix(),
        ...             tf.as_matrix(), atol=1e-12)
        True

        Multiple identity transformations can be generated at once:

        >>> tf = Tf.identity(2)
        >>> tf.as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]],
               [[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        >>> tf.single
        False
        >>> len(tf)
        2
        """
        if num is None:  # single
            return cls(np.eye(4, dtype=float), normalize=False, copy=False)
        else:
            return cls(np.eye(4, dtype=float)[None, :, :].repeat(num, axis=0),
                       normalize=False, copy=False)


    @cython.embedsignature(True)
    @classmethod
    def concatenate(cls, transformations):
        """
        Concatenate a sequence of `RigidTransformation` objects into a
        single object.

        Parameters
        ----------
        transformations : sequence of `RigidTransformation`
            If a single `RigidTransformation` instance is passed in, a copy of
            it is returned.

        Returns
        -------
        transformation : `RigidTransformation` instance
            The concatenated transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> tf1 = Tf.from_translation([1, 0, 0])
        >>> tf2 = Tf.from_translation([[2, 0, 0], [3, 0, 0]])
        >>> Tf.concatenate([tf1, tf2]).translation
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

        ..

            [R | t]
            [0 | 1]

        where ``R`` is a 3x3 orthonormal rotation matrix and ``t`` is a 3x1
        translation vector ``[tx, ty, tz]``.

        Returns
        -------
        matrix : numpy.ndarray, shape (4, 4) or (N, 4, 4)
            A single transformation matrix or a stack of transformation
            matrices.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        A transformation matrix is a 4x4 matrix formed from a 3x3 rotation
        matrix and a 3x1 translation vector:

        >>> t = np.array([2, 3, 4])
        >>> r = R.from_matrix([[0, 0, 1],
        ...                    [1, 0, 0],
        ...                    [0, 1, 0]])
        >>> tf = Tf.from_components(t, r)
        >>> tf.as_matrix()
        array([[ 0., 0., 1., 2.],
               [ 1., 0., 0., 3.],
               [ 0., 1., 0., 4.],
               [ 0., 0., 0., 1.]])

        >>> Tf.identity(2).as_matrix()
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
    def as_components(self):
        """Return the translation and rotation components of the transformation,
        where the rotation is applied first, followed by the translation.

        4x4 rigid transformation matrices are of the form:

        ..

            [R | t]
            [0 | 1]

        Where ``R`` is a 3x3 orthonormal rotation matrix and ``t`` is a 3x1
        translation vector ``[tx, ty, tz]``. This function returns the rotation
        corresponding to the this rotation matrix
        ``r = Rotation.from_matrix(R)`` and the translation vector ``t``.

        Take a transformation ``tf`` and a vector ``v``. When applying the
        transformation to the vector, the result is the same as if the
        transformation was applied to the vector in the following way:
        ``tf.apply(v) == translation + rotation.apply(v)``

        Returns
        -------
        translation : numpy.ndarray, shape (N, 3) or (3,)
            The translation of the transformation.
        rotation : `Rotation` instance
            The rotation of the transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Recover the rotation and translation from a transformation:

        >>> t = np.array([2, 3, 4])
        >>> r = R.from_matrix([[0, 0, 1],
        ...                    [1, 0, 0],
        ...                    [0, 1, 0]])
        >>> tf = Tf.from_components(t, r)
        >>> tf_t, tf_r = tf.as_components()
        >>> tf_t
        array([2., 3., 4.])
        >>> tf_r.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])

        The transformation applied to a vector is equivalent to the rotation
        applied to the vector followed by the translation:

        >>> r.apply([1, 0, 0])
        array([0., 1., 0.])
        >>> t + r.apply([1, 0, 0])
        array([2., 4., 4.])
        >>> tf.apply([1, 0, 0])
        array([2., 4., 4.])
        """
        return self.translation, self.rotation

    @cython.embedsignature(True)
    def as_exp_coords(self):
        """Return the exponential coordinates of the transformation.

        This implements the logarithmic map that converts SE(3) to
        6-dimensional real vectors. The first three components of the vector
        encode the rotation and the last three components encode the
        translation.

        Returns
        -------
        exp_coords : numpy.ndarray, shape (N, 6) or (6,)
            A single exponential coordinate vector or a stack of exponential
            coordinate vectors. The first three components define the
            rotation and the last three components define the translation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> import numpy as np

        Get exponential coordinates of the identity matrix:

        >>> Tf.identity().as_exp_coords()
        array([0., 0., 0., 0., 0., 0.])
        """
        exp_coords = np.empty((self._matrix.shape[0], 6), dtype=float)
        rotations = Rotation.from_matrix(self._matrix[:, :3, :3])
        rot_vec = rotations.as_rotvec()
        exp_coords[:, :3] = rot_vec
        exp_coords[:, 3:] = np.einsum('ijk,ik->ij',
                                      _compute_se3_log_translation_transform(rot_vec),
                                      self._matrix[:, :3, 3])

        if self._single:
            return exp_coords[0]
        else:
            return exp_coords

    @cython.embedsignature(True)
    def as_dual_quat(self, *, scalar_first=False):
        """Return the dual quaternion representation of the transformation.

        Unit dual quaternions encode orientation in a real unit quaternion
        and translation in a dual quaternion. There is a double cover, i.e.,
        the unit dual quaternions q and -q represent the same transformation.

        Parameters
        ----------
        scalar_first : bool, optional
            Whether the scalar component goes first or last in the two
            individual quaternions that represent the real and the dual part.
            Default is False, i.e. the scalar-last order is used.

        Returns
        -------
        dual_quat : numpy.ndarray, shape (N, 8) or (8,)
            A single unit dual quaternion vector or a stack of unit dual
            quaternion vectors. The real part is stored in the first four
            components and the dual part in the last four components.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> import numpy as np

        Get identity dual quaternion (we use scalar-last by default):

        >>> Tf.identity().as_dual_quat()
        array([0., 0., 0., 1., 0., 0., 0., 0.])

        When we want to use the scalar-first convention, we use the argument:

        >>> Tf.identity().as_dual_quat(scalar_first=True)
        array([1., 0., 0., 0., 0., 0., 0., 0.])
        """
        rotations = Rotation.from_matrix(self._matrix[:, :3, :3])
        real_parts = rotations.as_quat()

        pure_translation_quats = np.empty((len(self._matrix), 4), dtype=float)
        pure_translation_quats[:, :3] = self._matrix[:, :3, 3]
        pure_translation_quats[:, 3] = 0.0

        dual_parts = 0.5 * np.asarray(
            compose_quat(pure_translation_quats, real_parts))

        if scalar_first:
            real_parts = np.roll(real_parts, 1, axis=1)
            dual_parts = np.roll(dual_parts, 1, axis=1)

        dual_quats = np.hstack((real_parts, dual_parts))

        if self._single:
            return dual_quats[0]
        else:
            return dual_quats

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
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> tf = Tf.identity(3)
        >>> len(tf)
        3

        >>> tf = Tf.from_translation([1, 0, 0])
        >>> len(tf)  # doctest: +IGNORE_EXCEPTION_DETAIL
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
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> t = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]  # 3 translations
        >>> tf = Tf.from_translation(t)

        A single index returns a single transformation:

        >>> tf[0].as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])

        A slice returns a stack of transformations:

        >>> tf[1:3].translation
        array([[1., 0., 0.],
               [2., 0., 0.]])

        An index array returns a stack of transformations:

        >>> tf[[0, 2]].translation
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

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> t = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]  # 3 translations
        >>> tf = Tf.from_translation(t)

        Set a single transformation:

        >>> tf[0] = Tf.from_translation([9, 9, 9])
        >>> tf.translation
        array([[9., 9., 9.],
               [1., 0., 0.],
               [2., 0., 0.]])
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
        a vector ``v`` is equivalent to
        ``p.translation + p.rotation.apply(q.translation)
        + (p.rotation * q.rotation).apply(v)``.

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
              transformation in the other object with the single
              transformation.
            - Both ``p`` and ``q`` contain ``N`` transformations. In this case
              each transformation ``p[i]`` is composed with the corresponding
              transformation ``q[i]`` and `output` contains ``N``
              transformations.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Compose two transformations:

        >>> tf1 = Tf.from_translation([1, 0, 0])
        >>> tf2 = Tf.from_translation([0, 1, 0])
        >>> tf = tf1 * tf2
        >>> tf.translation
        array([1., 1., 0.])
        >>> tf.single
        True

        When applied to a vector, the composition of two transformations is
        applied in right-to-left order.

        >>> t1, r1 = [1, 2, 3], R.from_euler('z', 60, degrees=True)
        >>> t2, r2 = [0, 1, 0], R.from_euler('x', 30, degrees=True)
        >>> tf1 = Tf.from_components(t1, r1)
        >>> tf2 = Tf.from_components(t2, r2)
        >>> tf = tf1 * tf2
        >>> tf.apply([1, 0, 0])
        array([0.6339746, 3.3660254, 3.       ])
        >>> tf1.apply(tf2.apply([1, 0, 0]))
        array([0.6339746, 3.3660254, 3.       ])

        When at least one of the transformations is not single, the result is a
        stack of transformations.

        >>> tf1 = Tf.from_translation([1, 0, 0])
        >>> tf2 = Tf.from_translation([[0, 2, 0], [0, 0, 3]])
        >>> tf = tf1 * tf2
        >>> tf.translation
        array([[1., 2., 0.],
               [1., 0., 3.]])
        >>> tf.single
        False
        >>> len(tf)
        2
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
        itself `n` times. In other words,
        ``p ** -abs(n) == p.inv() ** abs(n)``.

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
        There are three notable cases: if ``n == 1`` then a copy of the
        original transformation is returned, if ``n == 0`` then the identity
        transformation is returned, and if ``n == -1`` then the inverse
        transformation is returned.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> import numpy as np

        A power of 2 returns the transformation composed with itself:

        >>> tf = Tf.from_translation([1, 2, 3])
        >>> (tf ** 2).translation
        array([2., 4., 6.])
        >>> (tf ** 2).as_matrix()
        array([[1., 0., 0., 2.],
               [0., 1., 0., 4.],
               [0., 0., 1., 6.],
               [0., 0., 0., 1.]])

        A negative power returns the inverse of the transformation raised to
        the absolute value of `n`:

        >>> (tf ** -2).translation
        array([-2., -4., -6.])
        >>> np.allclose((tf ** -2).as_matrix(), (tf.inv() ** 2).as_matrix(),
        ...             atol=1e-12)
        True

        A power of 0 returns the identity transformation:

        >>> (tf ** 0).as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])

        A power of 1 returns a copy of the original transformation:

        >>> (tf ** 1).as_matrix()
        array([[1., 0., 0., 1.],
               [0., 1., 0., 2.],
               [0., 0., 1., 3.],
               [0., 0., 0., 1.]])
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

            Tf = self.__class__(matrix, copy=False)
            if n < 0:
                return Tf.inv()
            else:
                return Tf

    @cython.embedsignature(True)
    def inv(self):
        """Invert this transformation.

        Composition of a transformation with its inverse results in an identity
        transformation.

        A rigid transformation is a composition of a rotation and a
        translation, where the rotation is applied first, followed by the
        translation. So the inverse transformation is equivalent to the inverse
        translation followed by the inverse rotation.

        Returns
        -------
        `RigidTransformation` instance
            The inverse of this transformation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        A transformation composed with its inverse results in an identity
        transformation:

        >>> rng = np.random.default_rng(seed=123)
        >>> t = rng.random(3)
        >>> r = R.random(rng=rng)
        >>> tf = Tf.from_components(t, r)
        >>> tf.as_matrix()
        array([[-0.45431291,  0.67276178, -0.58394466,  0.68235186],
               [-0.23272031,  0.54310598,  0.80676958,  0.05382102],
               [ 0.85990758,  0.50242162, -0.09017473,  0.22035987],
               [ 0.        ,  0.        ,  0.        ,  1.        ]])

        >>> (tf.inv() * tf).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])

        The inverse rigid transformation is the same as the inverse translation
        followed by the inverse rotation:

        >>> t, r = tf.as_components()
        >>> r_inv = r.inv()  # inverse rotation
        >>> t_inv = -t  # inverse translation
        >>> tf_r_inv = Tf.from_rotation(r_inv)
        >>> tf_t_inv = Tf.from_translation(t_inv)
        >>> np.allclose((tf_r_inv * tf_t_inv).as_matrix(),
        ...             tf.inv().as_matrix(),
        ...             atol=1e-12)
        True
        >>> (tf_r_inv * tf_t_inv * tf).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        """
        r_inv = np.swapaxes(self._matrix[:, :3, :3], 1, 2)  # Transpose to invert
        # This einsum performs element-wise matrix multiplication
        t_inv = -np.einsum('ijk,ik->ij', r_inv, self._matrix[:, :3, 3])

        matrix = _create_transformation_matrix(t_inv, r_inv, self._single)
        return self.__class__(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    def apply(self, vector, inverse=False):
        """Apply the transformation to a vector.

        If the original frame transforms to the final frame by this
        transformation, then its application to a vector can be seen in two
        ways:

            - As a projection of vector components expressed in the final frame
              to the original frame.
            - As the physical transformation of a vector being glued to the
              original frame as it transforms. In this case the vector
              components are expressed in the original frame before and after
              the transformation.

        In terms of rotation matrices and translation vectors, this application
        is the same as
        ``self.translation + self.rotation.as_matrix() @ vector``.

        Parameters
        ----------
        vector : array_like, shape (N, 3) or (3,)
            A single vector or a stack of vectors.
        inverse : bool, optional
            If True, the inverse of the transformation is applied to the
            vector.

        Returns
        -------
        transformed_vector : numpy.ndarray, shape (N, 3) or (3,)
            The transformed vector(s).
            Shape depends on the following cases:

                - If object contains a single transformation (as opposed to a
                  stack with a single transformation) and a single vector is
                  specified with shape ``(3,)``, then `transformed_vector` has
                  shape ``(3,)``.
                - In all other cases, `transformed_vector` has shape
                  ``(N, 3)``, where ``N`` is either the number of
                  transformations or vectors.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Apply a single transformation to a vector. Here the transformation is
        just a translation, so the result is the vector added to the
        translation vector.

        >>> t = np.array([1, 2, 3])
        >>> tf = Tf.from_translation(t)
        >>> t + np.array([1, 0, 0])
        array([2, 2, 3])
        >>> tf.apply([1, 0, 0])
        array([2., 2., 3.])

        Apply a single transformation to a stack of vectors:

        >>> tf.apply([[1, 0, 0], [0, 1, 0]])
        array([[2., 2., 3.],
               [1., 3., 3.]])

        Apply the inverse of a transformation to a vector, so the result is the
        negative of the translation vector added to the vector.

        >>> -t + np.array([1, 0, 0])
        array([0, -2, -3])
        >>> tf.apply([1, 0, 0], inverse=True)
        array([0., -2., -3.])

        For transformations which are not just pure translations, applying it
        to a vector is the same as applying the rotation component to the
        vector and then adding the translation component.

        >>> r = R.from_euler('z', 60, degrees=True)
        >>> tf = Tf.from_components(t, r)
        >>> t + r.apply([1, 0, 0])
        array([1.5,       2.8660254, 3.       ])
        >>> tf.apply([1, 0, 0])
        array([1.5,       2.8660254, 3.       ])

        When applying the inverse of a transformation, the result is the
        negative of the translation vector added to the vector, and then
        rotated by the inverse rotation.

        >>> r.inv().apply(-t + np.array([1, 0, 0]))
        array([-1.73205081, -1.        , -3.        ])
        >>> tf.apply([1, 0, 0], inverse=True)
        array([-1.73205081, -1.        , -3.        ])
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
        translated. This property returns the rotation part of the
        transformation.

        Returns
        -------
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        The rotation component is extracted from the transformation:

        >>> t = np.array([1, 0, 0])
        >>> r = R.random(3)
        >>> tf = Tf.from_components(t, r)
        >>> np.allclose(tf.rotation.as_matrix(), r.as_matrix())
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
        >>> from scipy.spatial.transform import RigidTransformation as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        The translation component is extracted from the transformation:

        >>> t = np.array([[1, 0, 0], [2, 0, 0], [3, 0, 0]])
        >>> r = R.random()
        >>> tf = Tf.from_components(t, r)
        >>> np.allclose(tf.translation, t)
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
