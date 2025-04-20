# cython: cpow=True

import numpy as np
from ._rotation_cy import compose_quat, from_matrix as quat_from_matrix, as_matrix as quat_as_matrix

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
def from_rotation(quat: double[:, :]) -> double[:, :, :]:
    rotmat = quat_as_matrix(quat)
    num_transforms = len(rotmat)
    matrix = np.zeros((num_transforms, 4, 4), dtype=float)
    matrix[:, :3, :3] = rotmat
    matrix[:, 3, 3] = 1
    return matrix


@cython.embedsignature(True)
def from_translation(translation) -> double[:, :, :]:
    translation = np.asarray(translation, dtype=float)

    if translation.ndim not in [1, 2] or translation.shape[-1] != 3:
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
    return matrix


@cython.embedsignature(True)
def from_components(translation, quat) -> double[:, :, :]:
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
