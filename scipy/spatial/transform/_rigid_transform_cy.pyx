# cython: cpow=True

import numpy as np
from ._rotation_cy import from_matrix as from_rot_matrix
from ._rotation_cy import _from_matrix_orthogonal as from_rot_matrix_orthogonal
from ._rotation_cy import inv as rot_inv
from ._rotation_cy import from_quat, from_rotvec
from ._rotation_cy import as_matrix, as_quat, as_rotvec, compose_quat
from ._rotation_cy import mean as rot_mean

cimport numpy as np
cimport cython

np.import_array()


@cython.embedsignature(True)
@cython.boundscheck(False)
def from_matrix(double[:, :, :] matrix, bint normalize=True, bint copy=True):
    mat = np.asarray(matrix, dtype=float)

    if mat.shape[-1] != 4 or mat.shape[-2] != 4:
        raise ValueError("Expected `matrix` to have shape (4, 4), or (N, 4, 4), "
                         f"got {mat.shape}.")

    if normalize or copy:
        mat = mat.copy()

    # Rigid transforms have the following matrix representation:
    # [R | t]
    # [0 | 1]
    # where R is a 3x3 orthonormal rotation matrix and t is a 3x1 translation
    # vector. The last row is always [0, 0, 0, 1] exactly

    # Check the last row. Exact match is required, as this last row should
    # not accumulate floating point error during composition.
    last_row_ok = np.all(mat[:, 3, :] == np.array([0, 0, 0, 1]), axis=1)
    if np.any(~last_row_ok):
        ind = np.where(~last_row_ok)[0][0]
        raise ValueError(f"Expected last row of transformation matrix {ind} to be "
                            f"exactly [0, 0, 0, 1], got {mat[ind, 3, :]}.")

    # The quat_from_matrix() method orthogonalizes the rotation
    # component of the transformation matrix. While this does have some
    # overhead in converting a rotation matrix to a quaternion and back, it
    # allows for skipping singular value decomposition for near-orthogonal
    # matrices, which is a computationally expensive operation.
    if normalize:
        mat[:, :3, :3] = as_matrix(from_rot_matrix(mat[:, :3, :3]))

    return mat


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def from_rotation(quat):
    single = quat.ndim == 1
    if quat.ndim == 1:
        quat = quat[None, :]
    # We don't need to normalize the quaternion here. The backend is only supposed to
    # be called from quaternions resulting from Rotation.as_quat(), which are already
    # guaranteed to be normalized.
    rotmat = as_matrix(from_quat(quat, normalize=False, copy=False))
    num_transforms = len(rotmat)
    matrix = np.zeros((num_transforms, 4, 4), dtype=float)
    matrix[:, :3, :3] = rotmat
    matrix[:, 3, 3] = 1
    if single:
        return matrix[0]
    return matrix


@cython.embedsignature(True)
@cython.boundscheck(False)
def from_translation(translation):
    translation = np.asarray(translation, dtype=float)

    if translation.ndim not in [1, 2] or translation.shape[-1] != 3:
        raise ValueError(
            f"Expected `translation` to have shape (..., 3), got {translation.shape}."
        )

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
@cython.boundscheck(False)
@cython.wraparound(False)
def from_exp_coords(exp_coords):
    exp_coords = np.asarray(exp_coords, dtype=float)
    single = exp_coords.ndim == 1

    exp_coords = np.atleast_2d(exp_coords)

    rot_vec = exp_coords[:, :3]
    rot_matrix = as_matrix(from_rotvec(rot_vec))
    translations = np.einsum('ijk,ik->ij',
                                _compute_se3_exp_translation_transform(rot_vec),
                                exp_coords[:, 3:])
    matrix = _create_transformation_matrix(translations, rot_matrix, False)
    if single:
        return matrix[0]
    return matrix


@cython.embedsignature(True)
@cython.boundscheck(False)
def from_dual_quat(dual_quat, *, bint scalar_first=False):
    dual_quat = np.asarray(dual_quat, dtype=float)

    if dual_quat.ndim not in [1, 2] or dual_quat.shape[-1] != 8:
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

    # _normalize_dual_quaternion() guarantees that the quaternion is already
    # normalized, so we don't need to normalize it again.
    translation = 2.0 * np.asarray(
        compose_quat(dual_part, rot_inv(real_part)))[:, :3]
    matrix = _create_transformation_matrix(translation, as_matrix(real_part), single)

    return matrix


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def as_exp_coords(double[:, :, :] matrix):
    exp_coords = np.empty((matrix.shape[0], 6), dtype=float)
    rot_vec = as_rotvec(from_rot_matrix_orthogonal(matrix[:, :3, :3]))
    exp_coords[:, :3] = rot_vec
    exp_coords[:, 3:] = np.einsum('ijk,ik->ij',
                                  _compute_se3_log_translation_transform(rot_vec),
                                  matrix[:, :3, 3])
    return exp_coords


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def as_dual_quat(double[:, :, :] matrix, *, bint scalar_first=False):
    real_parts = as_quat(from_rot_matrix_orthogonal(matrix[:, :3, :3]))

    pure_translation_quats = np.empty((len(matrix), 4), dtype=float)
    pure_translation_quats[:, :3] = matrix[:, :3, 3]
    pure_translation_quats[:, 3] = 0.0

    dual_parts = 0.5 * np.asarray(
        compose_quat(pure_translation_quats, real_parts))

    if scalar_first:
        real_parts = np.roll(real_parts, 1, axis=1)
        dual_parts = np.roll(dual_parts, 1, axis=1)

    dual_quats = np.hstack((real_parts, dual_parts))
    return dual_quats


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def compose_transforms(double[:, :, :] tf_matrix, double[:, :, :] other_tf_matrix):
    len_self = len(tf_matrix)
    len_other = len(other_tf_matrix)
    if not(len_self == 1 or len_other == 1 or len_self == len_other):
        raise ValueError("Expected equal number of transforms in both or a "
                            "single transform in either object, got "
                            f"{len_self} transforms in first and {len_other}"
                            "transforms in second object.")
    return np.matmul(tf_matrix, other_tf_matrix)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def inv(double[:, :, :] matrix):
    r_inv = np.swapaxes(matrix[:, :3, :3], 1, 2)  # Transpose to invert
    # This einsum performs element-wise matrix multiplication
    t_inv = -np.einsum('ijk,ik->ij', r_inv, matrix[:, :3, 3])
    matrix = _create_transformation_matrix(t_inv, r_inv, False)
    return np.asarray(matrix)


@cython.embedsignature(True)
@cython.boundscheck(False)
def apply(double[:, :, :] matrix, vector, bint inverse=False):
    vector = np.asarray(vector, dtype=float)
    if vector.ndim not in [1, 2] or vector.shape[-1] != 3:
        raise ValueError("Expected vector to have shape (N, 3), or (3,), "
                            f"got {vector.shape}.")

    if vector.ndim == 1:
        vector = vector[None, :]

    vector = np.hstack([vector, np.ones((vector.shape[0], 1))])

    if inverse:
        m = inv(matrix)
    else:
        m = matrix

    # This einsum performs matrix multiplication of each of the (4, 4)
    # matrices in `m` with the (4,) vectors in `vector`, with proper
    # broadcasting for different dimensions of `m` and `vector`.
    res = np.einsum('ijk,ik->ij', m, vector)[:, :3]
    return res


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pow(double[:, :, :] matrix, float n):
    # Exact short-cuts
    if n == 0:
        return np.tile(np.eye(4, dtype=float)[None, ...], (matrix.shape[0], 1, 1))
    elif n == -1:
        return inv(matrix)
    elif n == 1:
        return matrix
    return from_exp_coords(as_exp_coords(matrix) * n)


@cython.embedsignature(True)
@cython.boundscheck(False)
def mean(double[:, :, :] matrix, weights=None, axis=None):
    if matrix.shape[0] == 0:
        raise ValueError("Mean of an empty transform set is undefined.")

    # The Cython path assumes matrix is Nx4x4, so axis has to be None, 0, -1, (0,), (-1,),
    # or (). The code path is unchanged for any of the options except (), where we
    # immediately return the matrix
    if axis == ():
        return matrix

    if axis is None:
        axis = (0,)
    if isinstance(axis, int):
        axis = (axis,)
    if not isinstance(axis, tuple):  # Must be tuple by now
        raise ValueError("`axis` must be None, int, or tuple of ints.")
    if min(axis) < -1 or max(axis) > 0:
        raise ValueError(
            f"axis {axis} is out of bounds for transform with shape "
            f"{np.asarray(matrix).shape[:-2]}."
        )
    # Axis must be 0 for the cython backend. Everything else should have raised an
    # error during validation.
    axis = 0

    quat = as_quat(from_rot_matrix_orthogonal(matrix[:, :3, :3]))
    t = np.asarray(matrix[:, :3, 3])

    if weights is None:
        weights = np.ones(quat.shape[0])
    else:
        weights = np.asarray(weights)
        if np.any(weights < 0):
            raise ValueError("`weights` must be non-negative.")
        if weights.ndim != 1:
            raise ValueError(f"Expected `weights` to be 1 dimensional, got "
                             f"{weights.shape}.")
        if weights.shape[0] != matrix.shape[0]:
            raise ValueError("Expected `weights` to have number of values equal to "
                             f"number of transforms, got {weights.shape[0]} values and "
                             f"{matrix.shape[0]} transforms.")

    quat_mean = rot_mean(quat, weights=weights, axis=axis)
    t_mean = np.average(t, axis=axis, weights=weights)
    r_mean = as_matrix(np.asarray([quat_mean]))
    return _create_transformation_matrix(t_mean, r_mean, single=True)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def setitem(double[:, :, :] matrix, indexer, value):
    arr = np.asarray(matrix)
    arr[indexer] = value
    return arr


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


cdef _compute_se3_log_translation_transform(rot_vec):
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


@cython.embedsignature(True)
@cython.wraparound(False)
def normalize_dual_quaternion(double[:, :] dual_quat):
    """Normalize dual quaternion."""
    dq = np.asarray(dual_quat, dtype=float)
    real, dual = _normalize_dual_quaternion(
        dq[..., :4], dq[..., 4:])
    return np.concatenate((real, dual), axis=-1)


cdef _normalize_dual_quaternion(np.ndarray[double, ndim=2] real_part, np.ndarray[double, ndim=2] dual_part):
    """Ensure that the dual quaternion has unit norm.

    The norm is a dual number and must be 1 + 0 * epsilon, which means that
    the real quaternion must have unit norm and the dual quaternion must be
    orthogonal to the real quaternion.
    """
    real_part = np.copy(real_part)
    dual_part = np.copy(dual_part)

    real_norm = np.linalg.norm(real_part, axis=1)

    # special case: real quaternion is 0, we set it to identity
    zero_real_mask = real_norm == 0.0
    real_part[zero_real_mask, :4] = [0., 0., 0., 1.]
    real_norm[zero_real_mask] = 1.0

    # 1. ensure unit real quaternion
    real_part /= real_norm[:, np.newaxis]
    dual_part /= real_norm[:, np.newaxis]

    # 2. ensure orthogonality of real and dual quaternion
    dual_part -= np.sum(real_part * dual_part, axis=1)[:, np.newaxis] * real_part

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
