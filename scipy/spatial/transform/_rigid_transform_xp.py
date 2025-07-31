from types import EllipsisType

from scipy._lib._array_api import array_namespace, Array, is_lazy_array, xp_vector_norm
import scipy._lib.array_api_extra as xpx
from scipy.spatial.transform._rotation_xp import (
    as_matrix as quat_as_matrix,
    from_matrix as quat_from_matrix,
    from_rotvec as quat_from_rotvec,
    as_rotvec as quat_as_rotvec,
    compose_quat,
    from_quat,
    inv as quat_inv,
)
from scipy._lib.array_api_compat import device
from scipy._lib._util import broadcastable


def from_matrix(matrix: Array, normalize: bool = True, copy: bool = True) -> Array:
    xp = array_namespace(matrix)
    # Shape check should be done before calling this function
    if normalize or copy:
        matrix = xp.asarray(matrix, copy=True)

    last_row_ok = xp.all(
        matrix[..., 3, :] == xp.asarray([0, 0, 0, 1.0], device=device(matrix)), axis=-1
    )
    if is_lazy_array(matrix):
        matrix = xp.where(last_row_ok[..., None, None], matrix, xp.nan)
    elif xp.any(~last_row_ok):
        last_row_ok = xpx.atleast_nd(last_row_ok, ndim=1, xp=xp)
        matrix = xpx.atleast_nd(matrix, ndim=3, xp=xp)
        ind = int(xp.nonzero(~last_row_ok)[0][0])
        raise ValueError(
            f"Expected last row of transformation matrix {ind} to be "
            f"exactly [0, 0, 0, 1], got {matrix[ind, 3, xp.arange(4)]}."
        )
    # The quat_from_matrix() method orthogonalizes the rotation
    # component of the transformation matrix. While this does have some
    # overhead in converting a rotation matrix to a quaternion and back, it
    # allows for skipping singular value decomposition for near-orthogonal
    # matrices, which is a computationally expensive operation.
    if normalize:
        rotmat = quat_as_matrix(quat_from_matrix(matrix[..., :3, :3]))
        matrix = xpx.at(matrix)[..., :3, :3].set(rotmat)
    return matrix


def from_rotation(quat: Array) -> Array:
    xp = array_namespace(quat)
    rotmat = quat_as_matrix(quat)
    matrix = xp.zeros((*rotmat.shape[:-2], 4, 4), dtype=quat.dtype, device=device(quat))
    matrix = xpx.at(matrix)[..., :3, :3].set(rotmat)
    matrix = xpx.at(matrix)[..., 3, 3].set(1)
    return matrix


def from_translation(translation: Array) -> Array:
    xp = array_namespace(translation)

    if translation.ndim not in [1, 2] or translation.shape[-1] != 3:
        raise ValueError(
            "Expected `translation` to have shape (3,), or (N, 3), "
            f"got {translation.shape}."
        )
    _device = device(translation)
    eye = xp.eye(4, dtype=translation.dtype, device=_device)

    matrix = xpx.atleast_nd(eye, ndim=translation.ndim + 1, xp=xp)
    matrix = xp.zeros(
        (*translation.shape[:-1], 4, 4),
        dtype=translation.dtype,
        device=_device,
    )
    matrix = xpx.at(matrix)[...].set(xp.eye(4, dtype=translation.dtype, device=_device))
    matrix = xpx.at(matrix)[..., :3, 3].set(translation)
    return matrix


def from_components(translation: Array, quat: Array) -> Array:
    return compose_transforms(from_translation(translation), from_rotation(quat))


def from_exp_coords(exp_coords: Array) -> Array:
    rot_vec = exp_coords[..., :3]
    rot_matrix = quat_as_matrix(quat_from_rotvec(rot_vec))
    translation_transform = _compute_se3_exp_translation_transform(rot_vec)
    translations = (translation_transform @ exp_coords[..., 3:, None])[..., 0]
    return _create_transformation_matrix(translations, rot_matrix)


def from_dual_quat(dual_quat: Array, *, scalar_first: bool = False) -> Array:
    xp = array_namespace(dual_quat)

    if dual_quat.ndim not in [1, 2] or dual_quat.shape[-1] != 8:
        raise ValueError(
            "Expected `dual_quat` to have shape (8,), or (N, 8), "
            f"got {dual_quat.shape}."
        )

    real_part = dual_quat[..., :4]
    dual_part = dual_quat[..., 4:]
    if scalar_first:
        real_part = xp.roll(real_part, -1, axis=-1)
        dual_part = xp.roll(dual_part, -1, axis=-1)

    real_part, dual_part = _normalize_dual_quaternion(real_part, dual_part)

    rot_quat = from_quat(real_part)

    translation = 2.0 * compose_quat(dual_part, quat_inv(rot_quat))[..., :3]
    matrix = _create_transformation_matrix(translation, quat_as_matrix(rot_quat))

    return matrix


def as_exp_coords(matrix: Array) -> Array:
    xp = array_namespace(matrix)
    rot_vec = quat_as_rotvec(quat_from_matrix(matrix[..., :3, :3]))
    translation_transform = _compute_se3_log_translation_transform(rot_vec)
    translations = (translation_transform @ matrix[..., :3, 3][..., None])[..., 0]
    exp_coords = xp.concat([rot_vec, translations], axis=-1)
    return exp_coords


def as_dual_quat(matrix: Array, *, scalar_first: bool = False) -> Array:
    xp = array_namespace(matrix)
    real_parts = quat_from_matrix(matrix[..., :3, :3])

    pure_translation_quats = xp.empty(
        (*matrix.shape[:-2], 4), dtype=matrix.dtype, device=device(matrix)
    )
    pure_translation_quats = xpx.at(pure_translation_quats)[..., :3].set(
        matrix[..., :3, 3]
    )
    pure_translation_quats = xpx.at(pure_translation_quats)[..., 3].set(0.0)

    dual_parts = 0.5 * compose_quat(pure_translation_quats, real_parts)

    if scalar_first:
        real_parts = xp.roll(real_parts, 1, axis=-1)
        dual_parts = xp.roll(dual_parts, 1, axis=-1)

    dual_quats = xp.concat([real_parts, dual_parts], axis=-1)
    return dual_quats


def compose_transforms(tf_matrix: Array, other_tf_matrix: Array) -> Array:
    if not broadcastable(tf_matrix.shape, other_tf_matrix.shape):
        len_a = tf_matrix.shape[0] if tf_matrix.ndim == 3 else 1
        len_b = other_tf_matrix.shape[0] if other_tf_matrix.ndim == 3 else 1
        raise ValueError(
            "Expected equal number of transforms in both or a "
            "single transform in either object, got "
            f"{len_a} transforms in first and {len_b}"
            "transforms in second object."
        )
    return tf_matrix @ other_tf_matrix


def inv(matrix: Array) -> Array:
    xp = array_namespace(matrix)
    r_inv = xp.matrix_transpose(matrix[..., :3, :3])
    # Matrix multiplication of r_inv and translation vector
    t_inv = -(r_inv @ matrix[..., :3, 3][..., None])[..., 0]
    matrix = xp.zeros((*matrix.shape[:-2], 4, 4), dtype=matrix.dtype)
    matrix = xpx.at(matrix)[..., :3, :3].set(r_inv)
    matrix = xpx.at(matrix)[..., :3, 3].set(t_inv)
    matrix = xpx.at(matrix)[..., 3, 3].set(1)
    return matrix


def apply(matrix: Array, vector: Array, inverse: bool = False) -> Array:
    xp = array_namespace(matrix)
    if vector.ndim not in [1, 2] or vector.shape[-1] != 3:
        raise ValueError(
            f"Expected vector to have shape (N, 3), or (3,), got {vector.shape}."
        )
    vec = xp.empty((*vector.shape[:-1], 4), dtype=vector.dtype, device=device(vector))
    vec = xpx.at(vec)[..., :3].set(vector)
    vec = xpx.at(vec)[..., 3].set(1)
    vec = vec[..., None]

    if inverse:
        matrix = inv(matrix)

    # TODO: We raise a ValueError manually here because letting the function run its course would
    # raise heterogeneous error types and messages for different frameworks. However, the error only
    # mimics numpy's error message and does not provide the same amount of context.
    if not broadcastable(matrix.shape, vec.shape):
        raise ValueError("operands could not be broadcast together")

    # This einsum performs matrix multiplication of each of the (..., 4, 4) matrices in `matrix`
    # with the (..., 4) vectors in `vec`, with proper broadcasting for different dimensions.
    return (matrix @ vec)[..., :3, 0]


def pow(matrix: Array, n: float) -> Array:
    # Check if execution is eager. If so, we can branch and avoid computing all special cases
    xp = array_namespace(matrix)
    _device = device(matrix)
    if not is_lazy_array(matrix):
        if n == 0:
            identity = xp.eye(4, dtype=matrix.dtype, device=_device)
            identity = xpx.atleast_nd(identity, ndim=matrix.ndim, xp=xp)
            return xp.tile(identity, (*matrix.shape[:-2], 1, 1))
        elif n == -1:
            return inv(matrix)
        elif n == 1:
            return matrix
        return from_exp_coords(as_exp_coords(matrix) * n)
    # Lazy execution. We compute all special cases and the general case and combine them using
    # xp.where.
    result = from_exp_coords(as_exp_coords(matrix) * n)
    identity = xp.eye(4, dtype=matrix.dtype, device=_device)
    result = xp.where(n == 0, identity, result)
    result = xp.where(n == -1, inv(matrix), result)
    result = xp.where(n == 1, matrix, result)
    return result


def setitem(
    matrix: Array, indexer: int | slice | EllipsisType | None, value: Array
) -> Array:
    return xpx.at(matrix)[indexer].set(value)


def normalize_dual_quaternion(dual_quat: Array) -> Array:
    """Normalize dual quaternion."""
    xp = array_namespace(dual_quat)
    real, dual = _normalize_dual_quaternion(dual_quat[..., :4], dual_quat[..., 4:])
    return xp.concat((real, dual), axis=-1)


def _create_transformation_matrix(
    translations: Array, rotation_matrices: Array
) -> Array:
    if not translations.shape[:-1] == rotation_matrices.shape[:-2]:
        print(translations.shape, rotation_matrices.shape)
        raise ValueError(
            "The number of rotation matrices and translations must be the same."
        )
    xp = array_namespace(translations)
    matrix = xp.empty(
        (*translations.shape[:-1], 4, 4),
        dtype=translations.dtype,
        device=device(translations),
    )
    matrix = xpx.at(matrix)[..., :3, :3].set(rotation_matrices)
    matrix = xpx.at(matrix)[..., :3, 3].set(translations)
    matrix = xpx.at(matrix)[..., 3, :3].set(0)
    matrix = xpx.at(matrix)[..., 3, 3].set(1)
    return matrix


def _compute_se3_exp_translation_transform(rot_vec: Array) -> Array:
    """Compute the transformation matrix from the se3 translation part to SE3
    translation.

    The transformation matrix depends on the rotation vector.
    """
    xp = array_namespace(rot_vec)
    _device = device(rot_vec)
    _dtype = rot_vec.dtype
    angle = xp_vector_norm(rot_vec, axis=-1, keepdims=True, xp=xp)
    small_scale = angle < 1e-3

    k1_small = 0.5 - angle**2 / 24 + angle**4 / 720
    # Avoid division by zero for non-branching computations. The value will get discarded in the
    # xp.where selection.
    safe_angle = angle + xp.asarray(small_scale, dtype=_dtype, device=_device)
    k1 = (1.0 - xp.cos(angle)) / safe_angle**2
    k1 = xp.where(small_scale, k1_small, k1)

    k2_small = 1 / 6 - angle**2 / 120 + angle**4 / 5040
    # Again, avoid division by zero by adding one to all near-zero angles.
    safe_angle = angle + xp.asarray(small_scale, dtype=_dtype, device=_device)
    k2 = (angle - xp.sin(angle)) / safe_angle**3
    k2 = xp.where(small_scale, k2_small, k2)

    s = _create_skew_matrix(rot_vec)
    eye = xp.eye(3, dtype=_dtype, device=_device)

    return eye + k1[..., None] * s + k2[..., None] * s @ s


def _compute_se3_log_translation_transform(rot_vec: Array) -> Array:
    """Compute the transformation matrix from SE3 translation to the se3
    translation part.

    It is the inverse of `_compute_se3_exp_translation_transform` in a closed
    analytical form.
    """
    xp = array_namespace(rot_vec)
    _dtype = rot_vec.dtype
    _device = device(rot_vec)
    angle = xp_vector_norm(rot_vec, axis=-1, keepdims=True, xp=xp)
    mask = angle < 1e-3

    k_small = 1 / 12 + angle**2 / 720 + angle**4 / 30240
    safe_angle = angle + xp.asarray(mask, dtype=_dtype, device=_device)
    k = (1 - 0.5 * angle / xp.tan(0.5 * safe_angle)) / safe_angle**2
    k = xp.where(mask, k_small, k)

    s = _create_skew_matrix(rot_vec)

    return xp.eye(3, dtype=_dtype, device=_device) - 0.5 * s + k[..., None] * s @ s


def _create_skew_matrix(vec: Array) -> Array:
    """Create skew-symmetric (aka cross-product) matrix for stack of vectors."""
    xp = array_namespace(vec)
    result = xp.zeros((*vec.shape[:-1], 3, 3), dtype=vec.dtype, device=device(vec))
    result = xpx.at(result)[..., 0, 1].set(-vec[..., 2])
    result = xpx.at(result)[..., 0, 2].set(vec[..., 1])
    result = xpx.at(result)[..., 1, 0].set(vec[..., 2])
    result = xpx.at(result)[..., 1, 2].set(-vec[..., 0])
    result = xpx.at(result)[..., 2, 0].set(-vec[..., 1])
    result = xpx.at(result)[..., 2, 1].set(vec[..., 0])
    return result


def _normalize_dual_quaternion(
    real_part: Array, dual_part: Array
) -> tuple[Array, Array]:
    """Ensure that unit norm of the dual quaternion.

    The norm is a dual number and must be 1 + 0 * epsilon, which means that
    the real quaternion must have unit norm and the dual quaternion must be
    orthogonal to the real quaternion.
    """
    xp = array_namespace(real_part)

    real_norm = xp_vector_norm(real_part, axis=-1, keepdims=True, xp=xp)

    # special case: real quaternion is 0, we set it to identity
    zero_real_mask = real_norm == 0.0
    unit_quat = xp.asarray(
        [0.0, 0.0, 0.0, 1.0], dtype=real_part.dtype, device=device(real_part)
    )
    real_part = xp.where(zero_real_mask, unit_quat, real_part)
    real_norm = xp.where(zero_real_mask, 1.0, real_norm)

    # 1. ensure unit real quaternion
    real_part = real_part / real_norm
    dual_part = dual_part / real_norm

    # 2. ensure orthogonality of real and dual quaternion
    dual_part -= xp.sum(real_part * dual_part, axis=-1, keepdims=True) * real_part

    return real_part, dual_part
