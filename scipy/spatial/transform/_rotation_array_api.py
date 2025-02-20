from scipy._lib._array_api import array_namespace, Array
from scipy._lib.array_api_compat import device


def as_quat(
    quat: Array,
    *,
    normalize: bool = True,
    copy: bool = False,
    canonical: bool = False,
    scalar_first: bool = False,
) -> Array:
    # All operations are done non-branching to enable JIT-compilation
    if quat.shape[-1] != 4 or quat.shape[0] == 0:
        raise ValueError(f"Expected `quat` to have shape (..., 4), got {quat.shape}.")
    xp = array_namespace(quat)
    dtype = xp.float64 if quat.dtype == xp.float64 else xp.float32
    quat = xp.asarray(quat, dtype=dtype)
    _device = device(quat)
    scalar_first = xp.asarray(scalar_first, device=_device)
    normalize = xp.asarray(normalize, device=_device)
    copy = xp.asarray(copy, device=_device)
    canonical = xp.asarray(canonical, device=_device)

    quat = xp.where(scalar_first, xp.roll(quat, -1, axis=-1), quat)
    quat = xp.where(normalize | copy, xp.asarray(quat, copy=True), quat)
    quat_norm = xp.linalg.vector_norm(quat, axis=-1, keepdims=True)

    # TODO: This is a deviation from the behavior of the cython version because JIT code needs
    # to be strictly non-branching. If we have zero quaternions, we return NaNs
    # if xp.any(quat_norm == 0):
    #     raise ValueError("Found zero norm quaternions in `quat`.")
    quat_norm = xp.where(
        quat_norm == 0, xp.asarray([xp.nan], device=device(quat)), quat_norm
    )
    quat = xp.where(normalize, quat / quat_norm, quat)
    quat = xp.where(canonical, _quat_canonical(quat), quat)
    return quat


def from_quat(
    quat: Array,
    *,
    normalize: bool = True,
    copy: bool = False,
    scalar_first: bool = False,
) -> Array:
    # All operations are done non-branching to enable JIT-compilation
    if quat.shape[-1] != 4 or quat.shape[0] == 0:
        raise ValueError(f"Expected `quat` to have shape (..., 4), got {quat.shape}.")
    xp = array_namespace(quat)
    dtype = xp.float64 if quat.dtype == xp.float64 else xp.float32
    quat = xp.asarray(quat, dtype=dtype)
    _device = device(quat)
    scalar_first = xp.asarray(scalar_first, device=_device)
    normalize = xp.asarray(normalize, device=_device)
    copy = xp.asarray(copy, device=_device)

    quat = xp.where(scalar_first, xp.roll(quat, -1, axis=-1), quat)
    quat = xp.where(normalize | copy, xp.asarray(quat, copy=True), quat)
    quat_norm = xp.linalg.vector_norm(quat, axis=-1, keepdims=True)

    # TODO: This is a deviation from the behavior of the cython version because JIT code needs
    # to be strictly non-branching. If we have zero quaternions, we return NaNs
    # if xp.any(quat_norm == 0):
    #     raise ValueError("Found zero norm quaternions in `quat`.")
    quat_norm = xp.where(
        quat_norm == 0, xp.asarray([xp.nan], device=device(quat)), quat_norm
    )
    quat = xp.where(normalize, quat / quat_norm, quat)
    return quat


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
    matrix = xp.stack(matrix_elements, axis=-1).reshape((*quat.shape[:-1], 3, 3))
    return matrix


def apply(quat: Array, points: Array, inverse: bool = False) -> Array:
    xp = array_namespace(quat)
    subscripts = "...ji,...j->...i" if inverse else "...ij,...j->...i"
    return xp.einsum(subscripts, as_matrix(quat), points)
