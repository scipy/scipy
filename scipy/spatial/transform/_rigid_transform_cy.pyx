# cython: cpow=True

import numpy as np
from ._rotation_cy import from_rotvec, as_rotvec
from ._rotation_cy import from_matrix as quat_from_matrix, as_matrix as quat_as_matrix

cimport numpy as np
cimport cython

np.import_array()

@cython.embedsignature(True)
def from_matrix(matrix: double[:, :, :], normalize: bool = True, copy: bool = True) -> double[:, :, :]:
    matrix = np.asarray(matrix, dtype=float)

    if matrix.shape[-1] != 4 or matrix.shape[-2] != 4:
        raise ValueError("Expected `matrix` to have shape (4, 4), or (N, 4, 4), "
                         f"got {matrix.shape}.")

    if normalize or copy:
        matrix = matrix.copy()

    # Rigid transforms have the following matrix representation:
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

    # The quat_from_matrix() method orthogonalizes the rotation
    # component of the transformation matrix. While this does have some
    # overhead in converting a rotation matrix to a quaternion and back, it
    # allows for skipping singular value decomposition for near-orthogonal
    # matrices, which is a computationally expensive operation.
    if normalize:
        matrix[:, :3, :3] = quat_as_matrix(quat_from_matrix(matrix[:, :3, :3]))

    return matrix


@cython.embedsignature(True)
def from_rotation(quat: double[:] | double[:, :]) -> double[:, :] | double[:, :, :]:
    single = quat.ndim == 1
    if quat.ndim == 1:
        quat = quat[None, :]
    rotmat = quat_as_matrix(quat)
    num_transforms = len(rotmat)
    matrix = np.zeros((num_transforms, 4, 4), dtype=float)
    matrix[:, :3, :3] = rotmat
    matrix[:, 3, 3] = 1
    if single:
        return matrix[0]
    return matrix


@cython.embedsignature(True)
def from_translation(translation) -> double[:, :] | double[:, :, :]:
    translation = np.asarray(translation, dtype=float)

    if translation.ndim not in [1, 2] or translation.shape[-1] != 3:
        raise ValueError("Expected `translation` to have shape (3,), or (N, 3), "
                            f"got {translation.shape}.")

    # If a single translation vector is given, convert it to a 2D 1 x 3 matrix
    single = False
    if translation.ndim == 1:
        translation = translation[None, :]
        single = True

    num_translations = translation.shape[0]

    matrix = np.repeat(np.eye(4, dtype=float)[None, :, :], num_translations, axis=0)
    matrix[:, :3, 3] = translation

    if single:
        matrix = matrix[0]
    return matrix


@cython.embedsignature(True)
def from_components(translation: double[:] | double[:, :], quat: double[:] | double[:, :]) -> double[:, :] | double[:, :, :]:
    translation_matrix = from_translation(translation)
    single = translation.ndim == 1 and quat.ndim == 1
    if quat.ndim == 1:
        quat = quat[None, :]
    if translation_matrix.ndim == 2:
        translation_matrix = translation_matrix[None, :, :]
    tf = compose_transforms(translation_matrix, from_rotation(quat))
    if single:
        return tf[0]
    return tf


@cython.embedsignature(True)
def compose_transforms(tf_matrix: double[:, :, :], other_tf_matrix: double[:, :, :]) -> double[:, :, :]:
    len_self = len(tf_matrix)
    len_other = len(other_tf_matrix)
    if not(len_self == 1 or len_other == 1 or len_self == len_other):
        raise ValueError("Expected equal number of transforms in both or a "
                            "single transform in either object, got "
                            f"{len_self} transforms in first and {len_other}"
                            "transforms in second object.")
    return np.matmul(tf_matrix, other_tf_matrix)


@cython.embedsignature(True)
def from_exp_coords(exp_coords: double[:] | double[:, :]) -> double[:, :] | double[:, :, :]:
    exp_coords = np.asarray(exp_coords, dtype=float)
    single = exp_coords.ndim == 1

    exp_coords = np.atleast_2d(exp_coords)

    rot_vec = exp_coords[:, :3]
    rot_matrix = quat_as_matrix(from_rotvec(rot_vec))
    translations = np.einsum('ijk,ik->ij',
                                _compute_se3_exp_translation_transform(rot_vec),
                                exp_coords[:, 3:])
    matrix = _create_transformation_matrix(translations, rot_matrix, False)
    if single:
        return matrix[0]
    return matrix


@cython.embedsignature(True)
def as_exp_coords(matrix: double[:, :, :]) -> double[:, :]:
    exp_coords = np.empty((matrix.shape[0], 6), dtype=float)
    rot_vec = as_rotvec(quat_from_matrix(matrix[:, :3, :3]))
    exp_coords[:, :3] = rot_vec
    exp_coords[:, 3:] = np.einsum('ijk,ik->ij',
                                    _compute_se3_log_translation_transform(rot_vec),
                                    matrix[:, :3, 3])
    return exp_coords


@cython.embedsignature(True)
def inv(matrix: double[:, :, :]) -> double[:, :, :]:
    r_inv = np.swapaxes(matrix[:, :3, :3], 1, 2)  # Transpose to invert
    # This einsum performs element-wise matrix multiplication
    t_inv = -np.einsum('ijk,ik->ij', r_inv, matrix[:, :3, 3])
    matrix = _create_transformation_matrix(t_inv, r_inv, False)
    return matrix


def _create_transformation_matrix(
    translations: double[:] | double[:, :],
    rotation_matrices: double[:, :] | double[:, :, :],
    single: bool,
) -> double[:, :] | double[:, :, :]:
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
        raise ValueError(
            "The number of rotation matrices and translations must be the same."
        )

    matrix = np.empty((len(translations), 4, 4), dtype=float)
    matrix[:, :3, :3] = rotation_matrices
    matrix[:, :3, 3] = translations
    matrix[:, 3, :3] = 0.0
    matrix[:, 3, 3] = 1.0
    if single:
        return matrix[0]
    else:
        return matrix


cdef _compute_se3_exp_translation_transform(rot_vec):
    """Compute the transformation matrix from the se3 translation part to SE3
    translation.

    The transformation matrix depends on the rotation vector.
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


cdef _compute_se3_log_translation_transform(rot_vec: double[:, :]):
    """Compute the transformation matrix from SE3 translation to the se3
    translation part.

    It is the inverse of `_compute_se3_exp_translation_transform` in a closed
    analytical form.
    """
    angle = np.linalg.norm(rot_vec, axis=1)
    mask = angle < 1e-3

    k = np.empty(len(rot_vec))
    k[mask] = 1/12 + angle[mask]**2 / 720 + angle[mask]**4 / 30240
    k[~mask] = (1 - 0.5 * angle[~mask] / np.tan(0.5 * angle[~mask])) / angle[~mask]**2

    s = _create_skew_matrix(rot_vec)

    return np.identity(3) - 0.5 * s + k[:, None, None] * s @ s


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
