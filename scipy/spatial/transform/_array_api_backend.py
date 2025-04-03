"""Array API backend for the `Rotation` class.

This module provides generic, functional implementations of the `Rotation` class methods that work
with any Array API-compatible backend.
"""

# Parts of the implementation are adapted from the cython backend and
# https://github.com/jax-ml/jax/blob/d695aa4c63ffcebefce52794427c46bad576680c/jax/_src/scipy/spatial/transform.py.
import re

import numpy as np
from scipy._lib._array_api import (
    array_namespace,
    Array,
    ArrayLike,
    is_jax,
    xp_promote,
    xp_result_type,
    xp_vector_norm,
)
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


def from_matrix(matrix: Array) -> Array:
    xp = array_namespace(matrix)
    # Promote is not guaranteed to copy, but we want to ensure we don't modify the original array.
    matrix = xp.asarray(xp_promote(matrix, force_floating=True, xp=xp), copy=True)
    # DECISION: Left-handed case results in NaNs instead of raising an error. This is a deviation
    # from the cython implementation, which raises an error.
    mask = xp.linalg.det(matrix) <= 0
    mask = xp.broadcast_to(mask[..., None, None], mask.shape + (3, 3))
    matrix = xp.where(mask, xp.asarray(xp.nan, device=device(matrix)), matrix)

    gramians = matrix @ xp.matrix_transpose(matrix)
    # TODO: We need to orthogonalize only the non-orthogonal matrices, but jax does not support
    # non-concrete boolean indexing or any form of computation without statically known shapes. We
    # either have to branch depending on lazy/non-lazy frameworks or pay the performance penalty for
    # the SVD.
    is_orthogonal = xp.all(xpx.isclose(gramians, xp.eye(3), atol=1e-12, xp=xp))
    U, _, Vt = xp.linalg.svd(matrix)
    orthogonal_matrix = U @ Vt
    matrix = xp.where(is_orthogonal, matrix, orthogonal_matrix)

    matrix_trace = matrix[..., 0, 0] + matrix[..., 1, 1] + matrix[..., 2, 2]
    decision = xp.stack(
        [matrix[..., 0, 0], matrix[..., 1, 1], matrix[..., 2, 2], matrix_trace],
        axis=-1,
    )
    choice = xp.argmax(decision, axis=-1, keepdims=True)
    quat = xp.empty((*matrix.shape[:-2], 4), dtype=matrix.dtype)
    # While the Array API now does support advanced indexing, we still need to know the shape of all
    # arrays statically. This is not possible if we index based on the argmax, so we compute each
    # case and assemble the final result with `xp.where`.

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


def from_rotvec(rotvec: Array, degrees: bool = False) -> Array:
    xp = array_namespace(rotvec)
    rotvec = xp.asarray(xp_promote(rotvec, force_floating=True, xp=xp), copy=True)
    # TODO: Relax the shape check once we support proper broadcasting
    if rotvec.ndim not in [1, 2] or rotvec.shape[-1] != 3:
        raise ValueError(
            f"Expected `rot_vec` to have shape (3,) or (N, 3), got {rotvec.shape}"
        )
    degrees = xp.asarray(degrees, device=device(rotvec))
    rotvec = xp.where(degrees, _deg2rad(rotvec), rotvec)

    angle = xp_vector_norm(rotvec, axis=-1, keepdims=True)
    small_angle = angle <= 1e-3
    angle2 = angle**2
    small_scale = 0.5 - angle2 / 48 + angle2**2 / 3840
    # We need to handle the case where angle is 0 to avoid division by zero. We use the value of the
    # Taylor series approximation, but non-branching operations require that we still divide by the
    # angle. Since we do not use the result where the angle is close to 0, this is safe.
    div_angle = angle + xp.asarray(small_angle, dtype=angle.dtype)
    large_scale = xp.sin(angle / 2) / div_angle
    scale = xp.where(small_angle, small_scale, large_scale)
    quat = xp.concat([rotvec * scale, xp.cos(angle / 2)], axis=-1)
    return quat


def from_mrp(mrp: Array) -> Array:
    xp = array_namespace(mrp)
    mrp = xp.asarray(xp_promote(mrp, force_floating=True, xp=xp), copy=True)
    if mrp.ndim not in [1, 2] or mrp.shape[len(mrp.shape) - 1] != 3:
        raise ValueError(
            f"Expected `mrp` to have shape (3,) or (N, 3), got {mrp.shape}"
        )
    mrp2_plus_1 = xp.linalg.vecdot(mrp, mrp, axis=-1)[..., None] + 1
    q_no_norm = xp.concat([2 * mrp[..., :3], (2 - mrp2_plus_1)], axis=-1)
    quat = q_no_norm / mrp2_plus_1
    return quat


def from_euler(seq: str, angles: Array, degrees: bool = False) -> Array:
    xp = array_namespace(angles)
    num_axes = len(seq)
    if num_axes < 1 or num_axes > 3:
        raise ValueError(
            "Expected axis specification to be a non-empty "
            f"string of upto 3 characters, got {seq}"
        )

    intrinsic = re.match(r"^[XYZ]{1,3}$", seq) is not None
    extrinsic = re.match(r"^[xyz]{1,3}$", seq) is not None
    if not (intrinsic or extrinsic):
        raise ValueError(
            "Expected axes from `seq` to be from ['x', 'y', "
            f"'z'] or ['X', 'Y', 'Z'], got {seq}"
        )

    if any(seq[i] == seq[i + 1] for i in range(num_axes - 1)):
        raise ValueError(f"Expected consecutive axes to be different, got {seq}")

    angles = xp.asarray(xp_promote(angles, force_floating=True, xp=xp), copy=True)
    angles, is_single = _format_angles(angles, degrees, num_axes)
    axes = [_elementary_basis_index(x) for x in seq.lower()]
    q = _elementary_quat_compose(angles, axes, intrinsic)
    return q[0, ...] if is_single else q


def from_davenport(
    axes: Array, order: str, angles: Array | float, degrees: bool = False
) -> Array:
    xp = array_namespace(axes)
    if order in ["e", "extrinsic"]:  # Must be static, cannot be jitted
        extrinsic = True
    elif order in ["i", "intrinsic"]:
        extrinsic = False
    else:
        raise ValueError(
            "order should be 'e'/'extrinsic' for extrinsic sequences or 'i'/"
            f"'intrinsic' for intrinsic sequences, got {order}"
        )
    axes = xp.asarray(xp_promote(axes, force_floating=True, xp=xp), copy=True)
    # Angles could be a scalar, so we first need to convert it to an array before promoting
    angles = xp_promote(xp.asarray(angles, copy=True), force_floating=True, xp=xp)

    axes = xpx.atleast_nd(axes, ndim=2, xp=xp)
    if axes.shape[-1] != 3:
        raise ValueError("Axes must be vectors of length 3.")

    num_axes = axes.shape[0]
    if num_axes < 1 or num_axes > 3:
        raise ValueError("Expected up to 3 axes, got {}".format(num_axes))

    axes = axes / xp_vector_norm(axes, axis=-1, keepdims=True)

    # Check if axes are orthogonal. Checks are shape dependant and can therefore be jitted.
    axes_not_orthogonal = xp.zeros((*axes.shape[:-1], 1), dtype=xp.bool)
    if num_axes > 1:
        # Cannot be True yet, so we do not need to use xp.logical_or
        axes_not_orthogonal = xpx.at(axes_not_orthogonal)[..., 0].set(
            xp.abs(xp.vecdot(axes[..., 0, :], axes[..., 1, :])) > 1e-7
        )
    if num_axes > 2:
        axes_not_orthogonal = xpx.at(axes_not_orthogonal)[..., 0].set(
            xp.logical_or(
                xp.abs(xp.vecdot(axes[..., 1, :], axes[..., 2, :])) > 1e-7,
                axes_not_orthogonal[..., 0],
            )
        )

    axes = xp.where(axes_not_orthogonal, xp.asarray(xp.nan, device=device(axes)), axes)
    angles, single = _format_angles(angles, degrees, num_axes)

    dim_q = (4) if single else (angles.shape[0], 4)
    q = xp.zeros(dim_q, dtype=angles.dtype, device=device(angles))
    q = xpx.at(q)[..., 3].set(1)

    if not single:
        angles = angles[..., None]
    for i in range(num_axes):
        qi = from_rotvec(angles[:, i, ...] * axes[i, :])
        q = compose_quat(qi, q) if extrinsic else compose_quat(q, qi)
    return q


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


def as_rotvec(quat: Array, degrees: bool = False) -> Array:
    xp = array_namespace(quat)
    quat = _quat_canonical(quat)
    ax_norm = xp_vector_norm(quat[..., :3], axis=-1, keepdims=True)
    angle = 2 * xp.atan2(ax_norm, quat[..., 3][..., None])
    small_angle = angle <= 1e-3
    angle2 = angle**2
    small_scale = 2 + angle2 / 12 + 7 * angle2**2 / 2880
    # We need to handle the case where sin(angle/2) is 0 to avoid division by zero. We use the value
    # of the Taylor series approximation, but non-branching operations require that we still divide
    # by the sin. Since we do not use the result where the angle is close to 0, adding one to the
    # sin where we discard the result is safe.
    div_sin = xp.sin(angle / 2.0) + xp.asarray(small_angle, dtype=angle.dtype)
    large_scale = angle / div_sin
    scale = xp.where(small_angle, small_scale, large_scale)
    degrees = xp.asarray(degrees, device=device(quat))
    scale = xp.where(degrees, _rad2deg(scale), scale)
    rotvec = scale * quat[..., :3]
    return rotvec


def as_mrp(quat: Array) -> Array:
    xp = array_namespace(quat)
    one = xp.asarray(1.0, device=device(quat), dtype=quat.dtype)
    sign = xp.where(quat[..., 3, None] < 0, -one, one)
    denominator = 1.0 + sign * quat[..., 3, None]
    return sign * quat[..., :3] / denominator


def as_euler(quat: Array, seq: str, degrees: bool = False) -> Array:
    xp = array_namespace(quat)

    # Sanitize the sequence
    if len(seq) != 3:
        raise ValueError(f"Expected 3 axes, got {seq}.")
    intrinsic = re.match(r"^[XYZ]{1,3}$", seq) is not None
    extrinsic = re.match(r"^[xyz]{1,3}$", seq) is not None
    if not (intrinsic or extrinsic):
        raise ValueError(
            "Expected axes from `seq` to be from "
            "['x', 'y', 'z'] or ['X', 'Y', 'Z'], "
            "got {}".format(seq)
        )
    if any(seq[i] == seq[i + 1] for i in range(2)):
        raise ValueError(f"Expected consecutive axes to be different, got {seq}")

    _device = device(quat)
    axes = [_elementary_basis_index(x) for x in seq.lower()]
    axes = axes if extrinsic else axes[::-1]
    i, j, k = axes
    symmetric = i == k
    k = 3 - i - j if symmetric else k

    extrinsic = xp.asarray(extrinsic, device=_device)
    symmetric = xp.asarray(symmetric, device=_device)
    sign = xp.asarray(
        (i - j) * (j - k) * (k - i) // 2, dtype=quat.dtype, device=_device
    )
    # Permute quaternion elements
    a = xp.where(symmetric, quat[..., 3], quat[..., 3] - quat[..., j])
    b = xp.where(symmetric, quat[..., i], quat[..., i] + quat[..., k] * sign)
    c = xp.where(symmetric, quat[..., j], quat[..., j] + quat[..., 3])
    d = xp.where(symmetric, quat[..., k] * sign, quat[..., k] * sign - quat[..., i])

    angles = _get_angles(extrinsic, symmetric, sign, np.pi / 2, a, b, c, d)
    degrees = xp.asarray(degrees, device=_device)
    angles = xp.where(degrees, _rad2deg(angles), angles)
    return angles


def as_davenport(
    quat: Array, axes: ArrayLike, order: str, degrees: bool = False
) -> Array:
    xp = array_namespace(quat)
    axes = xp_promote(axes, force_floating=True, xp=xp)

    # Check argument validity
    if order in ["e", "extrinsic"]:
        extrinsic = True
    elif order in ["i", "intrinsic"]:
        extrinsic = False
    else:
        raise ValueError(
            "order should be 'e'/'extrinsic' for extrinsic "
            "sequences or 'i'/'intrinsic' for intrinsic "
            "sequences, got {}".format(order)
        )
    if axes.shape[0] != 3:
        raise ValueError(f"Expected 3 axes, got {axes.shape}.")
    if axes.shape[1] != 3:
        raise ValueError("Axes must be vectors of length 3.")

    # normalize axes
    axes = axes / xp_vector_norm(axes, axis=-1, keepdims=True)
    if xp.vecdot(axes[0, ...], axes[1, ...]) >= 1e-7:
        raise ValueError("Consecutive axes must be orthogonal.")
    if xp.vecdot(axes[1, ...], axes[2, ...]) >= 1e-7:
        raise ValueError("Consecutive axes must be orthogonal.")

    angles = _compute_davenport_from_quat(
        quat, axes[0, ...], axes[1, ...], axes[2, ...], extrinsic
    )
    if degrees:
        angles = _rad2deg(angles)
    return angles


def inv(quat: Array) -> Array:
    xp = array_namespace(quat)
    quat = xpx.at(quat)[..., :3].multiply(-1, copy=True, xp=xp)
    return quat


def magnitude(quat: Array) -> Array:
    xp = array_namespace(quat)
    sin_q = xp_vector_norm(quat[..., :3], axis=-1)
    cos_q = xp.abs(quat[..., 3])
    angles = 2 * xp.atan2(sin_q, cos_q)
    return angles


def approx_equal(
    quat: Array, other: Array, atol: float | None = None, degrees: bool = False
) -> Array:
    xp = array_namespace(quat)

    if not broadcastable(quat.shape, other.shape):
        raise ValueError(
            "Expected equal number of rotations in both "
            "or a single rotation in either object, "
            f"got {quat.shape[0]} rotations in first and {other.shape[0]} rotations in "
            "second object."
        )

    # DECISION: We cannot warn conditioned on the value of `degrees`. However, we should not need
    # to warn in the first place. If the user has set the degree flag and atol is None, the function
    # is still working as expected.
    atol = 1e-8 if atol is None else atol  # Value in radians
    atol = xp.asarray(atol, device=device(quat))
    degrees = xp.asarray(degrees, device=device(quat))
    atol = xp.where(degrees, _deg2rad(atol), atol)

    quat_result = compose_quat(other, inv(quat))
    angles = magnitude(quat_result)
    return angles < atol


def mean(quat: Array, weights: Array | None = None) -> Array:
    xp = array_namespace(quat)
    if quat.shape[0] == 0:
        raise ValueError("Mean of an empty rotation set is undefined.")

    # Branching code is okay for checks that include meta info such as shapes and types
    if weights is None:
        weights = xp.ones(
            quat.shape[:-1], dtype=xp_result_type(quat, force_floating=True, xp=xp)
        )
        weights = xpx.atleast_nd(weights, ndim=1, xp=xp)
    else:
        weights = xp.asarray(xp_promote(weights, force_floating=True, xp=xp), copy=True)
    # TODO: Missing full broadcasting support. We should relax this and only check if broadcasting
    # is possible.
    if weights.ndim != 1:
        raise ValueError(
            "Expected `weights` to be 1 dimensional, got shape {}.".format(
                weights.shape
            )
        )
    n_rot = quat.shape[0] if quat.ndim > 1 else 1
    if weights.shape[0] != n_rot:
        raise ValueError(
            "Expected `weights` to have number of values "
            "equal to number of rotations, got "
            "{} values and {} rotations.".format(weights.shape[0], n_rot)
        )
    # DECISION: We cannot check for negative weights because jit code needs to be non-branching. We
    # return NaN instead
    weights = xp.where(weights < 0, xp.asarray(xp.nan, device=device(quat)), weights)

    # Make sure we can transpose quat
    quat = xpx.atleast_nd(quat, ndim=2, xp=xp)
    K = (weights * quat.T) @ quat
    _, v = xp.linalg.eigh(K)
    return v[..., -1]


def reduce(
    quat: Array,
    left: Array | None = None,
    right: Array | None = None,
) -> tuple[Array, Array | None, Array | None]:
    # DECISION: We cannot have variable number of return arguments for jit compiled functions. We
    # therefore always return the indices, and filter out later.
    # TOOD: Properly support broadcasting.
    xp = array_namespace(quat)
    quat = xpx.atleast_nd(quat, ndim=2, xp=xp)
    if left is None:
        left = xp.ones_like(quat)
    if right is None:
        right = xp.ones_like(quat)

    if left is None and right is None:
        reduced = quat
        return reduced, None, None
    elif right is None:
        right = xp.asarray([[0.0, 0.0, 0.0, 1.0]])
    elif left is None:
        left = xp.asarray([[0.0, 0.0, 0.0, 1.0]])

    # Levi-Civita tensor for triple product computations
    e = xp.zeros((3, 3, 3), dtype=xp_result_type(quat, force_floating=True, xp=xp))
    e = xpx.at(e)[0, 1, 2].set(1)
    e = xpx.at(e)[1, 2, 0].set(1)
    e = xpx.at(e)[2, 0, 1].set(1)
    e = xpx.at(e)[0, 2, 1].set(-1)
    e = xpx.at(e)[2, 1, 0].set(-1)
    e = xpx.at(e)[1, 0, 2].set(-1)

    # We want to calculate the real components of q = l * p * r. It can
    # be shown that:
    #     qs = ls * ps * rs - ls * dot(pv, rv) - ps * dot(lv, rv)
    #          - rs * dot(lv, pv) - dot(cross(lv, pv), rv)
    # where ls and lv denote the scalar and vector components of l.

    p = quat
    ps, pv = _split_rotation(p, xp)
    ls, lv = _split_rotation(left, xp)
    rs, rv = _split_rotation(right, xp)

    # Compute each term without einsum (not accessible in the Array API)
    # First term: ls * ps * rs
    term1 = ls * ps * rs
    # Second term: ls * dot(pv, rv)
    term2 = ls * xp.sum(pv * rv, axis=-1)
    # Third term: ps * dot(lv, rv)
    term3 = ps * xp.sum(lv * rv, axis=-1)
    # Fourth term: rs * dot(lv, pv)
    term4 = rs * xp.sum(lv * pv, axis=-1)
    # Fifth term: dot(cross(lv, pv), rv)
    lv_cross_pv = xp.linalg.cross(lv, pv)
    term5 = xp.sum(lv_cross_pv * rv, axis=-1)

    qs = xp.abs(term1 - term2 - term3 - term4 - term5)

    # Find best indices from scalar components
    max_ind = xp.argmax(xp.reshape(qs, (qs.shape[0], -1)), axis=1)
    left_best = max_ind // rv.shape[0]
    right_best = max_ind % rv.shape[0]
    # Array API limitation: Integer index arrays are only allowed with integer indices
    # TODO: Can we somehow avoid this? Equivalent to x[left_best[0], :]
    all_idx = xp.arange(left.shape[-1])
    left = left[left_best[0], all_idx]
    right = right[right_best[0], all_idx]

    # Reduce the rotation using the best indices
    reduced = compose_quat(left, compose_quat(p, right))

    if left is None:
        left_best = None
    if right is None:
        right_best = None
    return reduced, left_best, right_best


def apply(quat: Array, points: Array, inverse: bool = False) -> Array:
    xp = array_namespace(quat)
    mat = as_matrix(quat)
    inverse = xp.asarray(inverse, device=device(quat))
    # We do not have access to einsum. To avoid broadcasting issues, we add a singleton dimension
    # to the points array and remove it after the operation.
    # TODO: We currently evaluate both branches of the where statement. For eager execution models,
    # this may significantly slow down the function. We should check that compilers can optimize
    # the where statement (e.g. in jax) and check if we can have an eager version that only
    # evaluates the branch that is needed.
    points = points[..., None]
    if not broadcastable(mat.shape, points.shape):
        raise ValueError(
            "Expected equal numbers of rotations and vectors "
            ", or a single rotation, or a single vector, got "
            f"{mat.shape[0]} rotations and {points.shape[0]} vectors."
        )

    return xp.where(inverse, mat.mT @ points, mat @ points)[..., 0]


def setitem(quat: Array, value: Array, indexer) -> Array:
    quat = xpx.at(quat)[indexer, ...].set(value)
    return quat


def align_vectors(
    a: Array, b: Array, weights: Array | None = None, return_sensitivity: bool = False
) -> tuple[Array, Array, Array]:
    xp = array_namespace(a)
    # Check input vectors
    a_original = xp_promote(a, force_floating=True, xp=xp)
    b_original = xp_promote(b, force_floating=True, xp=xp)
    # TODO: This function does not support broadcasting yet.
    a = xpx.atleast_nd(a_original, ndim=2, xp=xp)
    b = xpx.atleast_nd(b_original, ndim=2, xp=xp)
    if a.shape[-1] != 3:
        raise ValueError(
            "Expected input `a` to have shape (3,) or (N, 3), got {}".format(
                a_original.shape
            )
        )
    if b.shape[-1] != 3:
        raise ValueError(
            "Expected input `b` to have shape (3,) or (N, 3), got {}".format(
                b_original.shape
            )
        )
    if a.shape != b.shape:
        raise ValueError(
            "Expected inputs `a` and `b` to have same shapes"
            ", got {} and {} respectively.".format(a_original.shape, b_original.shape)
        )
    N = a.shape[0]

    # Check weights
    if weights is None:
        weights = xp.ones(N, device=device(a), dtype=a.dtype)
    else:
        weights = xp.asarray(weights, device=device(a), dtype=a.dtype)
        if weights.ndim != 1:
            raise ValueError(
                "Expected `weights` to be 1 dimensional, got shape {}.".format(
                    weights.shape
                )
            )
        if N > 1 and (weights.shape[0] != N):
            raise ValueError(
                "Expected `weights` to have number of values "
                "equal to number of input vectors, got "
                "{} values and {} vectors.".format(weights.shape[0], N)
            )
        # DECISION: We cannot check for negative weights because jit code needs to be non-branching.
        # We return NaN instead
        weights = xp.where(weights < 0, xp.asarray(xp.nan, device=device(a)), weights)

    # For the special case of a single vector pair, we use the infinite
    # weight code path
    weight_is_inf = xp.asarray([True]) if N == 1 else weights == xp.inf
    # DECISION: We cannot error out on multiple infinite weights. We return NaN instead.
    n_inf = xp.sum(xp.astype(weight_is_inf, a.dtype))
    weights = xp.where(n_inf > 1, xp.asarray(xp.nan, device=device(a)), weights)

    inf_branch = xp.any(weight_is_inf, axis=-1)
    # DECISION: We cannot compute both branches for all frameworks. There are two main reasons:
    # 1. Computing both for eager execution models is expensive.
    # 2. Some operations will fail when running the unused branch because of numerical and
    # algorithmical issues. Numpy e.g. will raise an exception when trying to compute the svd of a
    # matrix with infinite weights. To prevent this, we only compute the branch that is needed. Jax
    # jit however requires us to take the full compute graph. Therefore, we use xp.where for jax and
    # a branching version for non-jax frameworks. Might apply for other lazy frameworks as well.
    #
    # Note that we could also solve this by exploiting the externals of xpx.apply_where. However,
    # we'd have to rely on the implementation details of apply_where, which is something we should
    # avoid in the long run.
    if is_jax(xp):
        q_opt, rssd, sensitivity = _align_vectors(a, b, weights)
        q_opt_inf, rssd_inf, sensitivity_inf = _align_vectors_fixed(a, b, weights)
        q_opt = xp.where(~inf_branch, q_opt, q_opt_inf)
        rssd = xp.where(~inf_branch, rssd, rssd_inf)
        sensitivity = xp.where(~inf_branch, sensitivity, sensitivity_inf)
    else:
        if xp.any(inf_branch):
            q_opt, rssd, sensitivity = _align_vectors_fixed(a, b, weights)
        else:
            q_opt, rssd, sensitivity = _align_vectors(a, b, weights)
    return q_opt, rssd, sensitivity


def _align_vectors(a: Array, b: Array, weights: Array) -> tuple[Array, Array, Array]:
    xp = array_namespace(a)
    B = (weights[:, None] * a).mT @ b
    u, s, vh = xp.linalg.svd(B)

    # Correct improper rotation if necessary (as in Kabsch algorithm)
    neg_det = xp.linalg.det(u @ vh) < 0
    s = xpx.at(s)[..., -1].set(xp.where(neg_det, -s[..., -1], s[..., -1]))
    u = xpx.at(u)[..., :, -1].set(xp.where(neg_det, -u[..., :, -1], u[..., :, -1]))

    C = u @ vh

    # DECISION: We cannot branch on the condition because jit code needs to be non-branching. Hence,
    # we omit the check for uniqueness (s[1] + s[2] < 1e-16 * s[0])
    ssd = xp.sum(weights * xp.sum(b**2 + a**2, axis=1)) - 2 * xp.sum(s)
    rssd = xp.sqrt(xp.maximum(ssd, xp.zeros(1)))

    # TODO: We currently need to always compute the sensitivity matrix because jit code needs to be
    # non-branching. We should check if compilers can optimize the where statement (e.g. in jax)
    # and check if we can have an eager version that only evaluates the branch that is needed.
    # See xpx.apply_where, issue: https://github.com/data-apis/array-api-extra/pull/141
    zeta = (s[..., 0] + s[..., 1]) * (s[..., 1] + s[..., 2]) * (s[..., 2] + s[..., 0])
    kappa = s[..., 0] * s[..., 1] + s[..., 1] * s[..., 2] + s[..., 2] * s[..., 0]
    sensitivity = xp.mean(weights) / zeta * (kappa * xp.eye(3) + B @ B.mT)
    q_opt = from_matrix(C)
    return q_opt, rssd, sensitivity


def _align_vectors_fixed(
    a: Array, b: Array, weights: Array
) -> tuple[Array, Array, Array]:
    xp = array_namespace(a)
    N = a.shape[0]
    weight_is_inf = xp.asarray([True]) if N == 1 else weights == xp.inf
    # We cannot use boolean masks for indexing because of jax. For the same reason, we also cannot
    # use dynamic slices. As a workaround, we roll the array so that the infinitely weighted vector
    # is at index 0. We then use static slices to get the primary and secondary vectors.
    #
    # Note that argmax fulfils a secondary purpose here:
    # When we trace this function with jax, this function might get executed even if weight_is_inf
    # does not have a single valid entry. Argmax will still give us a valid index which allows us to
    # proceed with the function (the results are discarded anyways), whereas boolean indexing would
    # result in invalid, zero-shaped arrays.

    inf_idx = xp.argmax(xp.astype(weight_is_inf, xp.uint8))
    # Bug: torch.argmax returns a tensor, but does not support tensors as shifts in xp.roll. We
    # cannot convert to int because this raises a jax concretization error during jitting. This
    # will ideally be solved by an update of array-api-compat.
    # Tracking issue: https://github.com/data-apis/array-api/issues/914
    if not is_jax(xp):
        inf_idx = int(inf_idx)
    a_sorted = xp.roll(a, shift=-inf_idx, axis=0)
    b_sorted = xp.roll(b, shift=-inf_idx, axis=0)
    weights_sorted = xp.roll(weights, shift=-inf_idx, axis=0)

    a_pri = a_sorted[0, ...][None, ...]  # Effectively [[0], ...]
    b_pri = b_sorted[0, ...][None, ...]
    a_pri_norm = xp_vector_norm(a_pri, axis=-1, keepdims=True)
    b_pri_norm = xp_vector_norm(b_pri, axis=-1, keepdims=True)

    # We cannot error out on zero length vectors. We set the norm to NaN to avoid division by
    # zero and mark the result as invalid.
    a_pri_norm = xp.where(
        a_pri_norm == 0, xp.asarray(xp.nan, device=device(a)), a_pri_norm
    )
    b_pri_norm = xp.where(
        b_pri_norm == 0, xp.asarray(xp.nan, device=device(a)), b_pri_norm
    )

    a_pri = a_pri / a_pri_norm
    b_pri = b_pri / b_pri_norm

    # We first find the minimum angle rotation between the primary
    # vectors.
    cross = xp.linalg.cross(b_pri[..., 0, :], a_pri[..., 0, :])
    cross_norm = xp_vector_norm(cross, axis=-1)
    theta = xp.atan2(cross_norm, xp.sum(a_pri[..., 0, :] * b_pri[..., 0, :], axis=-1))
    tolerance = 1e-3  # tolerance for small angle approximation (rad)
    q_flip = xp.asarray([0.0, 0.0, 0.0, 1.0])

    # Near pi radians, the Taylor series approximation of x/sin(x) diverges, so for numerical
    # stability we flip pi and then rotate back by the small angle pi - theta
    flip = xp.asarray(np.pi - theta < tolerance)
    # For antiparallel vectors, cross = [0, 0, 0] so we need to manually set an arbitrary
    # orthogonal axis of rotation
    i = xp.argmin(xp.abs(a_pri[..., 0, :]))
    # TODO: Array API does not support fancy indexing __setitem__. The code is equivalent to doing
    # the following:
    # r_components = xp.zeros(3)
    # r_components = xpx.at(r_components)[i - 1].set(a_pri[0, i - 2])
    # r_components = xpx.at(r_components)[i - 2].set(-a_pri[0, i - 1])
    r_components = xp.stack(
        [
            xp.where(i == 1, a_pri[0, 2], xp.where(i == 2, -a_pri[0, 1], 0)),
            xp.where(i == 2, a_pri[0, 0], xp.where(i == 0, -a_pri[0, 2], 0)),
            xp.where(i == 0, a_pri[0, 1], xp.where(i == 1, -a_pri[0, 0], 0)),
        ]
    )
    r = xp.where(cross_norm == 0, r_components, cross)

    q_flip = xp.where(flip, from_rotvec(r / xp_vector_norm(r) * np.pi), q_flip)
    theta = xp.where(flip, np.pi - theta, theta)
    cross = xp.where(flip, -cross, cross)

    # Small angle Taylor series approximation for numerical stability
    theta2 = theta * theta
    small_scale = xp.abs(theta) < tolerance
    r_small_scale = cross * (1 + theta2 / 6 + theta2 * theta2 * 7 / 360)
    # We need to handle the case where theta is 0 to avoid division by zero. We use the value of the
    # Taylor series approximation, but non-branching operations require that we still divide by the
    # angle. Since we do not use the result where the angle is close to 0, this is safe.
    theta = theta + xp.asarray(small_scale, dtype=theta.dtype)
    r_large_scale = cross * theta / xp.sin(theta)
    r = xp.where(small_scale, r_small_scale, r_large_scale)
    q_pri = compose_quat(from_rotvec(r), q_flip)

    # Case 1): No secondary vectors, q_opt is q_pri -> Immediately done
    # Case 2): Secondary vectors exist
    # We cannot use boolean masks here because of jax
    a_sec = a_sorted[1:, ...]
    b_sec = b_sorted[1:, ...]
    weights_sec = weights_sorted[1:]

    # We apply the first rotation to the b vectors to align the
    # primary vectors, resulting in vectors c. The secondary
    # vectors must now be rotated about that primary vector to best
    # align c to a.
    c_sec = apply(q_pri, b_sec)

    # Calculate vector components for the angle calculation. The
    # angle phi to rotate a single 2D vector C to align to 2D
    # vector A in the xy plane can be found with the equation
    # phi = atan2(cross(C, A), dot(C, A))
    #     = atan2(|C|*|A|*sin(phi), |C|*|A|*cos(phi))
    # The below equations perform the same operation, but with the
    # 3D rotation restricted to the 2D plane normal to a_pri, where
    # the secondary vectors are projected into that plane. We then
    # find the composite angle from the weighted sum of the
    # axial components in that plane, minimizing the 2D alignment
    # problem.
    sin_term = xp.linalg.vecdot(xp.linalg.cross(c_sec, a_sec), a_pri)
    cos_term = xp.linalg.vecdot(c_sec, a_sec) - (
        xp.linalg.vecdot(c_sec, a_pri) * xp.linalg.vecdot(a_sec, a_pri)
    )
    phi = xp.atan2(xp.sum(weights_sec * sin_term), xp.sum(weights_sec * cos_term))
    q_sec = from_rotvec(phi * a_pri[0, ...])

    # Compose these to get the optimal rotation
    q_opt = xp.where(xp.asarray(N == 1), q_pri, compose_quat(q_sec, q_pri))

    # Calculated the root sum squared distance. We force the error to
    # be zero for the infinite weight vectors since they will align
    # exactly.
    weights_inf_zero = xp.asarray(weights, copy=True)

    multiple_vectors = xp.asarray(N > 1, device=device(weights))
    mask = xp.logical_or(multiple_vectors, weights[0] == xp.inf)
    mask = xp.logical_and(mask, weight_is_inf)
    # Skip non-infinite weight single vectors pairs, we used the
    # infinite weight code path but don't want to zero that weight
    weights_inf_zero = xpx.at(weights_inf_zero)[mask].set(0)
    a_est = apply(q_opt, b)
    rssd = xp.sqrt(xp.sum(weights_inf_zero @ (a - a_est) ** 2))

    mask = xp.any(xp.isnan(weights), axis=-1)
    q_opt = xp.where(mask, xp.asarray(xp.nan, device=device(q_opt)), q_opt)
    return q_opt, rssd, xp.asarray(xp.nan, device=device(q_opt))


def pow(quat: Array, n: float) -> Array:
    xp = array_namespace(quat)
    # general scaling of rotation angle
    result = from_rotvec(n * as_rotvec(quat))
    # Special cases 0 -> identity, -1 -> inv, 1 -> copy
    identity = xp.zeros(
        (*quat.shape[:-1], 4), dtype=result.dtype, device=device(result)
    )
    identity = xpx.at(identity)[..., 3].set(1)
    mask = xp.asarray(n == 0, device=device(quat))
    result = xp.where(mask, identity, result)
    mask = xp.asarray(n == -1, device=device(quat))
    result = xp.where(mask, inv(quat), result)
    mask = xp.asarray(n == 1, device=device(quat))
    result = xp.where(mask, quat, result)
    return result


def _normalize_quaternion(quat: Array) -> Array:
    xp = array_namespace(quat)
    quat_norm = xp_vector_norm(quat, axis=-1, keepdims=True)
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


def _format_angles(angles: Array, degrees: bool, num_axes: int) -> tuple[Array, bool]:
    xp = array_namespace(angles)
    angles = xp.where(xp.asarray(degrees), _deg2rad(angles), angles)

    is_single = False
    # Prepare angles to have shape (num_rot, num_axes)
    if num_axes == 1:
        if angles.ndim == 0:
            # (1, 1)
            angles = xpx.atleast_nd(angles, ndim=2, xp=xp)
            is_single = True
        elif angles.ndim == 1:
            # (N, 1)
            angles = angles[:, None]
        elif angles.ndim == 2 and angles.shape[-1] != 1:
            raise ValueError(
                "Expected `angles` parameter to have shape (N, 1), got {}.".format(
                    angles.shape
                )
            )
        elif angles.ndim > 2:
            raise ValueError(
                "Expected float, 1D array, or 2D array for "
                "parameter `angles` corresponding to `seq`, "
                "got shape {}.".format(angles.shape)
            )
    else:  # 2 or 3 axes
        if angles.ndim not in [1, 2] or angles.shape[-1] != num_axes:
            raise ValueError(
                "Expected `angles` to be at most "
                "2-dimensional with width equal to number "
                "of axes specified, got "
                "{} for shape".format(angles.shape)
            )

        if angles.ndim == 1:
            # (1, num_axes)
            angles = angles[None, :]
            is_single = True

    # By now angles should have shape (num_rot, num_axes)
    # sanity check
    if angles.ndim != 2 or angles.shape[-1] != num_axes:
        raise ValueError(
            "Expected angles to have shape (num_rotations, num_axes), got {}.".format(
                angles.shape
            )
        )

    return angles, is_single


def _compute_davenport_from_quat(
    quat: Array, n1: Array, n2: Array, n3: Array, extrinsic: bool
) -> Array:
    # The algorithm assumes extrinsic frame transformations. The algorithm
    # in the paper is formulated for rotation quaternions, which are stored
    # directly by Rotation.
    # Adapt the algorithm for our case by reversing both axis sequence and
    # angles for intrinsic rotations when needed
    xp = array_namespace(quat)
    mask = xp.asarray(extrinsic, device=device(quat))
    _n1 = xp.where(mask, n1, n3)
    _n3 = xp.where(mask, n3, n1)
    n1, n3 = _n1, _n3

    n_cross = xp.linalg.cross(n1, n2)
    lamb = xp.atan2(xp.vecdot(n3, n_cross), xp.vecdot(n3, n1))

    # alternative set of angles compatible with as_euler implementation
    mask = xp.asarray(lamb < 0, device=device(quat))
    n2 = xp.where(mask, -n2, n2)
    lamb = xp.where(mask, -lamb, lamb)
    n_cross = xp.where(mask, -n_cross, n_cross)
    correct_set = xp.where(mask, xp.asarray(True), xp.asarray(False))

    quat_lamb = xp.asarray([*(xp.sin(lamb / 2) * n2), xp.cos(lamb / 2)])

    q_trans = compose_quat(quat_lamb, quat)
    a = q_trans[..., 3]
    b = xp.linalg.vecdot(q_trans[..., :3], n1)
    c = xp.linalg.vecdot(q_trans[..., :3], n2)
    d = xp.linalg.vecdot(q_trans[..., :3], n_cross)

    angles = _get_angles(extrinsic, False, 1, lamb, a, b, c, d)
    angles = xpx.at(angles)[..., 1].set(
        xp.where(correct_set, -angles[..., 1], angles[..., 1])
    )

    return angles


def _elementary_quat_compose(angles: Array, axes: list[int], intrinsic: bool) -> Array:
    xp = array_namespace(angles)
    intrinsic = xp.asarray(intrinsic, device=device(angles))
    quat = _make_elementary_quat(axes[0], angles[..., 0])
    for i in range(1, len(axes)):
        ax_quat = _make_elementary_quat(axes[i], angles[..., i])
        quat = xp.where(
            intrinsic, compose_quat(quat, ax_quat), compose_quat(ax_quat, quat)
        )
    return quat


def _make_elementary_quat(axis: int, angle: Array) -> Array:
    xp = array_namespace(angle)
    quat = xp.zeros(
        (*angle.shape, 4), dtype=xp_result_type(angle, force_floating=True, xp=xp)
    )
    quat = xpx.at(quat)[..., 3].set(xp.cos(angle / 2.0))
    quat = xpx.at(quat)[..., axis].set(xp.sin(angle / 2.0))
    return quat


def _get_angles(
    extrinsic: bool,
    symmetric: bool,
    sign: int,
    lamb: float,
    a: Array,
    b: Array,
    c: Array,
    d: Array,
) -> Array:
    xp = array_namespace(a)
    eps = 1e-7
    half_sum = xp.atan2(b, a)
    half_diff = xp.atan2(d, c)
    angles = xp.zeros((*a.shape, 3), dtype=a.dtype)
    _device = device(a)

    angles = xpx.at(angles)[..., 1].set(2 * xp.atan2(xp.hypot(c, d), xp.hypot(a, b)))

    angle_first = 0 if extrinsic else 2
    angle_third = 2 if extrinsic else 0

    # Convert extrinsic and symmetric to arrays for use in xp.where
    extrinsic = xp.asarray(extrinsic, dtype=xp.bool, device=_device)
    symmetric = xp.asarray(symmetric, dtype=xp.bool, device=_device)

    # Check if the second angle is close to 0 or pi, causing a singularity.
    # - Case 0: Second angle is neither close to 0 nor pi.
    # - Case 1: Second angle is close to 0.
    # - Case 2: Second angle is close to pi.
    case1 = xp.abs(angles[..., 1]) <= eps
    case2 = xp.abs(angles[..., 1] - np.pi) <= eps
    case0 = ~(case1 | case2)

    one = xp.asarray(1, dtype=angles.dtype, device=_device)
    a0 = xp.where(case1, 2 * half_sum, 2 * half_diff * xp.where(extrinsic, -one, one))
    angles = xpx.at(angles)[..., 0].set(a0)

    a1 = xp.where(case0, half_sum - half_diff, angles[..., angle_first])
    angles = xpx.at(angles)[..., angle_first].set(a1)

    a3 = xp.where(case0, half_sum + half_diff, angles[..., angle_third])
    a3 = xp.where(symmetric, a3, a3 * sign)
    angles = xpx.at(angles)[..., angle_third].set(a3)

    a1 = xp.where(symmetric, angles[..., 1], angles[..., 1] - lamb)
    angles = xpx.at(angles)[..., 1].set(a1)

    angles = (angles + np.pi) % (2 * np.pi) - np.pi

    return angles


def compose_quat(p: Array, q: Array) -> Array:
    xp = array_namespace(p)
    cross = xp.linalg.cross(p[..., :3], q[..., :3])
    qx = p[..., 3] * q[..., 0] + q[..., 3] * p[..., 0] + cross[..., 0]
    qy = p[..., 3] * q[..., 1] + q[..., 3] * p[..., 1] + cross[..., 1]
    qz = p[..., 3] * q[..., 2] + q[..., 3] * p[..., 2] + cross[..., 2]
    qw = (
        p[..., 3] * q[..., 3]
        - p[..., 0] * q[..., 0]
        - p[..., 1] * q[..., 1]
        - p[..., 2] * q[..., 2]
    )
    quat = xp.stack([qx, qy, qz, qw], axis=-1)
    return quat


def broadcastable(shape_a: tuple[int, ...], shape_b: tuple[int, ...]) -> bool:
    """Check if two shapes are broadcastable."""
    return all(
        (m == n) or (m == 1) or (n == 1) for m, n in zip(shape_a[::-1], shape_b[::-1])
    )


def _split_rotation(q, xp):
    q = xpx.atleast_nd(q, ndim=2, xp=xp)
    return q[..., -1], q[..., :-1]


def _deg2rad(angles: Array) -> Array:
    return angles * np.pi / 180.0


def _rad2deg(angles: Array) -> Array:
    return angles * 180.0 / np.pi
