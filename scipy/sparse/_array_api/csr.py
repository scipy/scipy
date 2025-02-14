from scipy import sparse
import numpy as np
import numpy.typing as npt

def asarray(obj, *, dtype=None, copy=None, device=None):
    # This is a HUGE hack to get the array-api test
    # suite even running, much less doing anything sensible.
    if np.isscalar(obj):
         obj = np.array(obj).reshape(1, -1)
    if isinstance(obj, list):
        obj = np.array(obj).reshape(len(obj), -1)
    return sparse.csr_array(obj, dtype=dtype, copy=copy)

def empty(shape, *, dtype=None, device=None):
    return sparse.csr_array(shape, dtype=dtype)

def equal(x1, x2, /):
    return x1 == x2

def empty_like(x, *, dtype=None, device=None):
    return sparse.csr_array(x, dtype=dtype, copy=False)

def eye(n_rows, n_cols=None, /, *, k=0, dtype=None, device=None):
    return sparse.eye_array(n_rows, n_cols, k=k, dtype=dtype, format="csr")

def full(shape, fill_value, *, dtype=None, device=None):
    return sparse.csr_array(np.full(shape, fill_value, dtype=dtype))

def full_like(x, fill_value, *, dtype=None, device=None):
    return sparse.csr_array(np.full_like(x, fill_value, dtype=dtype))

def ones(shape, *, dtype=None, device=None):
    return sparse.csr_array(np.ones(shape, dtype=dtype))

def ones_like(x, *, dtype=None, device=None):
    return sparse.csr_array(np.ones_like(x, dtype=dtype))

def zeros(shape, *, dtype=None, device=None):
    if isinstance(shape, int):
        return sparse.csr_array((shape,1), dtype=dtype)
    return sparse.csr_array(shape, dtype=dtype)

def zeros_like(x, *, dtype=None, device=None):
    if dtype is None:
        dtype = x.dtype
    return sparse.csr_array(x.shape, dtype=dtype)

def arange(start, stop, step, *, dtype=None, device=None):
    return sparse.csr_array(np.arange(start, stop, step, dtype=dtype))

def linspace(start, stop, /, num, *, dtype=None, device=None, endpoint=True):
    return sparse.csr_array(
        np.linspace(start, stop, num, dtype=dtype, endpoint=endpoint)
    )

def any(x, /, *, axis = None, keepdims=False):
    if axis is None:
        return x.count_nonzero()
    return x.astype('bool').sum(axis=axis)

def less(x1, x2):
    return x1 < x2

def less_equal(x1, x2):
    return x1 <= x2

def greater(x1, x2):
    return x1 > x2

def greater_equal(x1, x2):
    return x1 >= x2

def floor(x, /):
    return sparse.csr_array((np.floor(x.data), x.indices, x.indptr), shape=x.shape)

def take(x, indices, /, *, axis: int | None = None):
    if axis is None:
        if x.ndim > 1:
            raise ValueError("axis must be specified when input is > 1d")
        return x[indices]
    idxer = [slice(None)] * x.ndim
    idxer[axis] = indices
    return x[tuple(idxer)]

def all(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False): # type: ignore
    return np.sum(x, axis=axis)

def isnan(x, /):
    return sparse.csr_array((np.isnan(x.data), x.indices, x.indptr), shape=x.shape)

def isinf(x, /):
    return sparse.csr_array((np.isinf(x.data), x.indices, x.indptr), shape=x.shape)

def reshape(x, /, shape: tuple[int, ...], *, copy: bool | None = None):
    raise NotImplementedError

def isfinite(x, /):
    return sparse.csr_array((np.isfinite(x.data), x.indices, x.indptr), shape=x.shape)

def mean(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.mean(axis=axis)

def min(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.min(axis=axis)

def max(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.max(axis=axis)

def sum(
        x,
        /,
        *,
        axis: int | tuple[int, ...] | None = None,
        dtype: npt.DTypeLike | None = None,
        keepdims: bool = False
    ):
    return x.sum(axis=axis)

def argmax(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.argmax(axis=axis)

def argmin(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.argmin(axis=axis)

def var(
    x,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    correction: int | float = 0.0,
    keepdims: bool = False
):
    divisor = 1
    if isinstance(axis, tuple):
        if len(axis) != 1:
            raise NotImplementedError(
                "Cannot calculate variance over more than one axis for sparse"
            )
        axis = axis[0]
        divisor = (x.shape[axis] - correction)
    return (1 / divisor) * (x ** 2).mean(axis=axis) - np.square(x.mean(axis=axis))

def std(
        x,
        /,
        *,
        axis: int | tuple[int, ...] | None = None,
        correction: int | float = 0.0,
        keepdims: bool = False
    ):
    return np.sqrt(var(x, axis=axis, correction=correction, keepdims=keepdims))

def unique_values(x, /):
    return np.append(np.unique(x.data), np.array([0], dtype=x.data.dtype))

def unique_counts(x, /):
    is_first_element_zero = x[0, 0] == 0
    values, counts = np.unique(x.data, return_counts=True)
    number_of_zeros = (x.shape[0] * x.shape[1]) - x.count_nonzero()
    if is_first_element_zero:
        values = np.append(np.array([0], dtype=x.data.dtype), values)
        counts = np.append(np.array([number_of_zeros]), values)
    else:
        values = np.append(values, np.array([0], dtype=x.data.dtype))
        counts = np.append(values, np.array([number_of_zeros]))
    return values, counts

def unique_all(x, /):
    raise NotImplementedError(
        'Construction of `inverse_indices` for sparse arrays is inefficient'
    )

def unique_inverse(x, /):
    raise NotImplementedError(
        'Construction of `inverse_indices` for sparse arrays is inefficient'
    )

def where(condition, x1, x2, /):
    data = x2.copy()
    data[condition] = x1
    return data

def logical_and(x1, x2, /):
    x1_bool = x1.astype('bool') 
    x2_bool = x2.astype('bool')
    return x1_bool * x2_bool

def logical_or(x1, x2, /):
    x1_bool = x1.astype('bool') 
    x2_bool = x2.astype('bool')
    return x1_bool + x2_bool

def logical_xor(x1, x2, /):
    x1_bool = x1.astype('bool') 
    x2_bool = x2.astype('bool')
    return x2_bool != x1_bool

def multiply(x1, x2, /):
    # TODO: Does scipy obey the type promotions listed for multiply?
    # https://data-apis.org/array-api/latest/API_specification/generated/array_api.multiply.html#multiply
    return x1.multiply(x2)

def subtract(x1, x2, /):
    # TODO: Does scipy obey type promotion rules in general?
    # https://data-apis.org/array-api/latest/API_specification/type_promotion.html#type-promotion
    return x1 - x2

def matmul(x1, x2, /):
    # TODO: In general, this will rely on 1D array functionality via broadcasting:
    # https://data-apis.org/array-api/latest/API_specification/generated/array_api.matmul.html#array_api.matmul
    return x1 @ x2

def not_equal(x1, x2, /):
    return x1 != x2

def binopt(x1, x2, op, /):
    x2_nonzero = np.column_stack(x2.nonzero())
    x1_nonzero = np.column_stack(x1.nonzero())
    non_zero_coords = np.moveaxis(
        np.unique(np.concatenate((x2_nonzero, x1_nonzero)), axis=0), 1, 0
    )
    x1_data = x1[non_zero_coords][0]
    x2_data = x2[non_zero_coords][0]
    op_data = np.ravel(op(x1_data, x2_data))
    return sparse.coo_array((op_data, (non_zero_coords[0], non_zero_coords[1]))).tocsr()

def bitwise_left_shift(x1, x2, /):
    shifted_data = np.left_shift(x1.data, x2)
    return sparse.csr_array((shifted_data, x1.indices, x1.indptr), shape=x1.shape)

def bitwise_right_shift(x1, x2, /):
    shifted_data = np.right_shift(x1.data, x2)
    return sparse.csr_array((shifted_data, x1.indices, x1.indptr), shape=x1.shape)

def bitwise_invert(x, /):
    inverted_data = np.bitwise_not(x.data)
    return sparse.csr_array((inverted_data, x.indices, x.indptr), shape=x.shape)

def bitwise_or(x1, x2, /):
    binopt(x1, x2, np.bitwise_or)

def bitwise_and(x1, x2, /):
    binopt(x1, x2, np.bitwise_and)
