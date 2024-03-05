from scipy import sparse
import numpy as np
import numpy.typing as npt

def asarray(obj, *, dtype=None, copy=None):
    # This is a HUGE hack to get the array-api test suite even running, much less doing anything sensible.
    if np.isscalar(obj):
         obj = np.array(obj).reshape(1, -1)
    if isinstance(obj, list):
        obj = np.array(obj).reshape(len(obj), -1)
    return sparse.csr_array(obj, dtype=dtype, copy=copy)

def empty(shape, *, dtype=None):
    return sparse.csr_array(shape, dtype=dtype)

def empty_like(x, *, dtype=None):
    return sparse.csr_array(x, dtype=dtype, copy=False)

def eye(n_rows, n_cols=None, /, *, k=0, dtype=None):
    return sparse.eye(n_rows, n_cols, k=k, dtype=dtype, format="csr")

def full(shape, fill_value, *, dtype=None):
    return sparse.csr_array(np.full(shape, fill_value, dtype=dtype))

def full_like(x, fill_value, *, dtype=None):
    return sparse.csr_array(np.full_like(x, fill_value, dtype=dtype))

def ones(shape, *, dtype=None):
    return sparse.csr_array(np.ones(shape, dtype=dtype))

def ones_like(x, *, dtype=None):
    return sparse.csr_array(np.ones_like(x, dtype=dtype))

def zeros(shape, *, dtype=None):
    if isinstance(shape, int):
        return sparse.csr_array((shape,1), dtype=dtype)
    return sparse.csr_array(shape, dtype=dtype)

def zeros_like(x, *, dtype=None):
    if dtype is None:
        dtype = x.dtype
    return sparse.csr_array(x.shape, dtype=dtype)

def floor(x):
    return sparse.csr_array((np.floor(x.data), x.indices, x.indptr), shape=x.shape)

def take(x, indices, /, *, axis: int | None = None):
    idxer = [slice(None)] * x.ndim
    idxer[axis] = indices
    return x[tuple(idxer)]

def all(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False): # type: ignore
    return np.sum(x, axis=axis)

def isnan(x):
    return sparse.csr_array((np.isnan(x.data), x.indices, x.indptr), shape=x.shape)

def isinf(x):
    return sparse.csr_array((np.isinf(x.data), x.indices, x.indptr), shape=x.shape)

def reshape(x, /, shape: tuple[int, ...], *, copy: bool | None = None):
    raise NotImplementedError

def isfinite(x):
    return sparse.csr_array((np.isfinite(x.data), x.indices, x.indptr), shape=x.shape)

def mean(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.mean(axis=axis)

def min(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.min(axis=axis)

def max(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.max(axis=axis)

def sum(x, /, *, axis: int | tuple[int, ...] | None = None, dtype: npt.DTypeLike | None = None, keepdims: bool = False):
    return x.sum(axis=axis)

def argmax(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.argmax(axis=axis)

def argmin(x, /, *, axis: int | tuple[int, ...] | None = None, keepdims: bool = False):
    return x.argmin(axis=axis)
