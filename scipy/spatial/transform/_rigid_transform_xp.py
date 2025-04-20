from scipy._lib._array_api import (
    array_namespace,
    Array,
    ArrayLike,
    is_lazy_array,
    xp_promote,
    xp_result_type,
    xp_vector_norm,
)
import scipy._lib.array_api_extra as xpx
from scipy.spatial.transform._rotation_xp import (
    as_matrix as quat_as_matrix,
    from_matrix as quat_from_matrix,
)
from scipy._lib.array_api_compat import device


def from_matrix(matrix: Array, normalize: bool = True, copy: bool = True) -> Array:
    xp = array_namespace(matrix)
    # Shape check should be done before calling this function
    if normalize or copy:
        matrix = xp.asarray(matrix, copy=True)

    last_row_ok = xp.all(matrix[..., 3, :] == xp.asarray([0, 0, 0, 1.0]))
    if is_lazy_array(matrix):
        matrix = xp.where(last_row_ok, matrix, xp.nan)
    elif xp.any(~last_row_ok):
        last_row_ok = xpx.atleast_nd(last_row_ok, ndim=1, xp=xp)
        matrix = xpx.atleast_nd(matrix, ndim=3, xp=xp)
        ind = xp.nonzero(~last_row_ok)[0][0]
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
    num_transforms = rotmat.shape[0]
    matrix = xp.zeros((num_transforms, 4, 4), dtype=quat.dtype)
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
    matrix = xpx.atleast_nd(
        xp.eye(4, dtype=translation.dtype), ndim=translation.ndim + 1, xp=xp
    )
    _device = device(translation)
    matrix = xp.zeros(
        (*translation.shape[:-1], 4, 4),
        dtype=translation.dtype,
        device=_device,
    )
    matrix = xpx.at(matrix)[...].set(xp.eye(4, dtype=translation.dtype, device=_device))
    matrix = xpx.at(matrix)[..., :3, 3].set(translation)
    return matrix
