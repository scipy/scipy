import re

import numpy as np
from scipy._lib._array_api import array_namespace, Array, ArrayLike
from scipy._lib.array_api_compat import device
import scipy._lib.array_api_extra as xpx


def from_quat(
    quat: Array,
    normalize: bool = True,
    copy: bool = True,
    *,
    scalar_first: bool = False,
) -> Array:
    xp = array_namespace(quat)
    _device = device(quat)
    normalize = xp.asarray(normalize, device=_device)
    copy = xp.asarray(copy, device=_device)
    scalar_first = xp.asarray(scalar_first, device=_device)
    quat = xp.where(scalar_first, xp.roll(quat, -1, axis=-1), quat)
    quat = xp.where(normalize | copy, xp.asarray(quat, copy=True), quat)
    quat = xp.where(normalize, _normalize_quaternion(quat), quat)
    return quat


def from_euler(seq: str, angles: Array, degrees: bool = False) -> Array:
    xp = array_namespace(angles)
    num_axes = len(seq)
    if num_axes < 1 or num_axes > 3:
        raise ValueError(
            "Expected axis specification to be a non-empty "
            "string of upto 3 characters, got {}".format(seq)
        )

    intrinsic = re.match(r"^[XYZ]{1,3}$", seq) is not None
    extrinsic = re.match(r"^[xyz]{1,3}$", seq) is not None
    if not (intrinsic or extrinsic):
        raise ValueError(
            "Expected axes from `seq` to be from ['x', 'y', "
            "'z'] or ['X', 'Y', 'Z'], got {}".format(seq)
        )

    if any(seq[i] == seq[i + 1] for i in range(num_axes - 1)):
        raise ValueError(
            "Expected consecutive axes to be different, got {}".format(seq)
        )

    dtype = xp.float32
    if angles.dtype != xp.float32:
        dtype = xp.result_type(xp.float32, xp.float64)
    angles = xp.asarray(angles, dtype=dtype)
    angles = xpx.atleast_nd(angles, ndim=1, xp=xp)
    axes = xp.asarray([_elementary_basis_index(x) for x in seq.lower()])
    return _elementary_quat_compose(angles, axes, intrinsic, degrees)


def from_matrix(matrix: Array) -> Array:
    xp = array_namespace(matrix)
    matrix = xp.asarray(matrix, copy=True, dtype=atleast_f32(matrix))
    # DECISION: Left-handed case results in NaNs instead of raising an error. This is a deviation
    # from the cython implementation, which raises an error.
    # TODO: Masking is only necessary because the Array API does not yet support advanced indexing
    # with arrays of indices. See comment further below.
    mask = xp.linalg.det(matrix) <= 0
    mask_shape = (1,) * (len(matrix.shape) - 2) + (3, 3)
    mask = xp.tile(mask[..., None, None], mask_shape)
    matrix = xpx.at(matrix)[mask].set(xp.nan)

    gramians = matrix @ xp.matrix_transpose(matrix)
    # TODO: Masking is only necessary because the Array API does not yet support advanced indexing
    # with arrays of indices. See comment further below.
    is_orthogonal = xpx.isclose(gramians, xp.eye(3), atol=1e-12, xp=xp)
    mask = ~xp.all(is_orthogonal, axis=(-2, -1))
    mask = xp.tile(mask[..., None, None], mask_shape)
    non_orthogonal_matrix = xp.reshape(matrix[mask], (-1, 3, 3))
    U, _, Vt = xp.linalg.svd(non_orthogonal_matrix)
    orthogonal_matrix = U @ Vt
    matrix_flat = xp.reshape(matrix, (-1,))
    mask_flat = xp.reshape(mask, (-1,))
    orthogonal_matrix_flat = xp.reshape(orthogonal_matrix, (-1,))
    matrix_flat = xpx.at(matrix_flat)[mask_flat].set(orthogonal_matrix_flat)
    matrix = xp.reshape(matrix_flat, matrix.shape)

    matrix_trace = matrix[..., 0, 0] + matrix[..., 1, 1] + matrix[..., 2, 2]
    decision = xp.stack(
        [matrix[..., 0, 0], matrix[..., 1, 1], matrix[..., 2, 2], matrix_trace],
        axis=-1,
    )
    choice = xp.argmax(decision, axis=-1, keepdims=True)
    quat = xp.empty((*matrix.shape[:-2], 4), dtype=matrix.dtype)
    # TODO: The Array API does not yet support advanced indexing with arrays of indices, so we
    # compute each case and assemble the final result with `xp.where`. Advanced indexing is
    # currently under development, see https://github.com/data-apis/array-api/issues/669.
    # As soon as this makes it into the spec, we can optimize this function.
    # https://github.com/data-apis/array-api/milestone/4

    # Case 0
    quat_0 = xp.stack(
        [
            1 - matrix_trace[...] + 2 * matrix[..., 0, 0],
            matrix[..., 1, 0] + matrix[..., 0, 1],
            matrix[..., 2, 0] + matrix[..., 0, 2],
            matrix[..., 2, 1] - matrix[..., 1, 2],
        ],
        axis=-1,
    )
    quat = xp.where((choice == 0), quat_0, quat)

    # Case 1
    quat_1 = xp.stack(
        [
            matrix[..., 1, 0] + matrix[..., 0, 1],
            1 - matrix_trace[...] + 2 * matrix[..., 1, 1],
            matrix[..., 2, 1] + matrix[..., 1, 2],
            matrix[..., 0, 2] - matrix[..., 2, 0],
        ],
        axis=-1,
    )
    quat = xp.where((choice == 1), quat_1, quat)

    # Case 2
    quat_2 = xp.stack(
        [
            matrix[..., 2, 0] + matrix[..., 0, 2],
            matrix[..., 2, 1] + matrix[..., 1, 2],
            1 - matrix_trace[...] + 2 * matrix[..., 2, 2],
            matrix[..., 1, 0] - matrix[..., 0, 1],
        ],
        axis=-1,
    )
    quat = xp.where((choice == 2), quat_2, quat)

    # Case 3
    quat_3 = xp.stack(
        [
            matrix[..., 2, 1] - matrix[..., 1, 2],
            matrix[..., 0, 2] - matrix[..., 2, 0],
            matrix[..., 1, 0] - matrix[..., 0, 1],
            1 + matrix_trace[...],
        ],
        axis=-1,
    )
    quat = xp.where((choice == 3), quat_3, quat)

    return _normalize_quaternion(quat)


def as_quat(
    quat: Array, canonical: bool = False, *, scalar_first: bool = False
) -> Array:
    xp = array_namespace(quat)
    _device = device(quat)
    scalar_first = xp.asarray(scalar_first, device=_device)
    canonical = xp.asarray(canonical, device=_device)
    quat = xp.where(canonical, _quat_canonical(quat), quat)
    quat = xp.where(scalar_first, xp.roll(quat, 1, axis=-1), quat)
    return quat


def as_matrix(quat: Array) -> Array:
    xp = array_namespace(quat)
    x = quat[..., 0]
    y = quat[..., 1]
    z = quat[..., 2]
    w = quat[..., 3]

    x2 = x * x
    y2 = y * y
    z2 = z * z
    w2 = w * w

    xy = x * y
    zw = z * w
    xz = x * z
    yw = y * w
    yz = y * z
    xw = x * w

    matrix_elements = [
        x2 - y2 - z2 + w2,
        2 * (xy - zw),
        2 * (xz + yw),
        2 * (xy + zw),
        -x2 + y2 - z2 + w2,
        2 * (yz - xw),
        2 * (xz - yw),
        2 * (yz + xw),
        -x2 - y2 + z2 + w2,
    ]
    matrix = xp.reshape(xp.stack(matrix_elements, axis=-1), (*quat.shape[:-1], 3, 3))
    return matrix


def apply(quat: Array, points: Array, inverse: bool = False) -> Array:
    xp = array_namespace(quat)
    subscripts = "...ji,...j->...i" if inverse else "...ij,...j->...i"
    return xp.einsum(subscripts, as_matrix(quat), points)


def inv(quat: Array) -> Array:
    xp = array_namespace(quat)
    quat = xpx.at(quat)[..., :3].multiply(-1, copy=True, xp=xp)
    return quat


def _normalize_quaternion(quat: Array) -> Array:
    xp = array_namespace(quat)
    quat_norm = xp.linalg.vector_norm(quat, axis=-1, keepdims=True)
    quat = xp.where(quat_norm == 0, xp.asarray([xp.nan], device=device(quat)), quat)
    return quat / quat_norm


def _quat_canonical(quat: Array) -> Array:
    xp = array_namespace(quat)
    mask = quat[..., 3] < 0
    zero_w = quat[..., 3] == 0
    mask = xp.logical_or(mask, zero_w & (quat[..., 0] < 0))
    zero_wx = xp.logical_or(zero_w, quat[..., 0] == 0)
    mask = xp.logical_or(mask, zero_wx & (quat[..., 1] < 0))
    zero_wxy = xp.logical_or(zero_wx, quat[..., 1] == 0)
    mask = xp.logical_or(mask, zero_wxy & (quat[..., 2] < 0))
    return xp.where(mask[..., None], -quat, quat)


def _elementary_basis_index(axis: str) -> int:
    if axis == "x":
        return 0
    elif axis == "y":
        return 1
    elif axis == "z":
        return 2
    raise ValueError(f"Expected axis to be from ['x', 'y', 'z'], got {axis}")


def _elementary_quat_compose(
    angles: Array, axes: Array, intrinsic: bool, degrees: bool
) -> Array:
    xp = array_namespace(angles)
    degrees = xp.asarray(degrees, device=device(angles))
    intrinsic = xp.asarray(intrinsic, device=device(angles))
    angles = xp.where(degrees, _deg2rad(angles), angles)
    quat0 = _make_elementary_quat(axes[0], angles[..., 0])
    quat1 = _make_elementary_quat(axes[1], angles[..., 1])
    quat2 = _make_elementary_quat(axes[2], angles[..., 2])
    quat = xp.where(intrinsic, compose_quat(quat0, quat1), compose_quat(quat1, quat0))
    quat = xp.where(intrinsic, compose_quat(quat, quat2), compose_quat(quat2, quat))
    return quat


def _make_elementary_quat(axis: int, angle: Array) -> Array:
    xp = array_namespace(angle)
    quat = xp.zeros((*angle.shape[1:], 4), dtype=atleast_f32(angle))
    quat = xpx.at(quat)[..., 3].set(xp.cos(angle / 2.0))
    quat = xpx.at(quat)[..., axis].set(xp.sin(angle / 2.0))
    return quat


def compose_quat(p: Array, q: Array) -> Array:
    xp = array_namespace(p)
    cross = xp.linalg.cross(p[..., :3], q[..., :3])
    return xp.asarray(
        [
            p[..., 3] * q[..., 0] + q[..., 3] * p[..., 0] + cross[..., 0],
            p[..., 3] * q[..., 1] + q[..., 3] * p[..., 1] + cross[..., 1],
            p[..., 3] * q[..., 2] + q[..., 3] * p[..., 2] + cross[..., 2],
            p[..., 3] * q[..., 3]
            - p[..., 0] * q[..., 0]
            - p[..., 1] * q[..., 1]
            - p[..., 2] * q[..., 2],
        ]
    )


def _deg2rad(angles: Array) -> Array:
    return angles * np.pi / 180.0


def atleast_f32(x: ArrayLike) -> type:
    xp = array_namespace(x)
    # In case it's an Array and it's float32, we do not promote
    if getattr(x, "dtype", None) == xp.float32:
        return xp.float32
    return xp.result_type(xp.float32, xp.float64)
