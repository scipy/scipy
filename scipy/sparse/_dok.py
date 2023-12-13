"""Dictionary Of Keys based matrix"""

__docformat__ = "restructuredtext en"

__all__ = ['dok_array', 'dok_matrix', 'isspmatrix_dok']

import itertools
import numpy as np

from ._matrix import spmatrix
from ._base import _spbase, sparray, issparse
from ._index import IndexMixin
from ._sputils import (
    isdense,
    getdtype,
    isshape,
    isintlike,
    isscalarlike,
    upcast,
    upcast_scalar,
    check_shape,
)


class _dok_base(_spbase, IndexMixin):
    _format = 'dok'

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        _spbase.__init__(self)
        self._dict = sd = {}

        self.dtype = getdtype(dtype, default=float)
        is_array = isinstance(self, sparray)
        if isinstance(arg1, tuple):
            if isshape(arg1, allow_1d=is_array):
                self._shape = check_shape(arg1, allow_1d=is_array)
        elif issparse(arg1):  # Sparse ctor
            if arg1.format == self.format:
                if copy:
                    arg1 = arg1.copy()
            else:
                arg1 = arg1.todok()

            if dtype is not None:
                arg1 = arg1.astype(dtype, copy=False)

            sd.update(arg1._dict)
            self._shape = check_shape(arg1.shape, allow_1d=is_array)
            self.dtype = arg1.dtype
        else:  # Dense ctor
            try:
                arg1 = np.asarray(arg1)
            except Exception as e:
                raise TypeError('Invalid input format.') from e

            if arg1.ndim > 2:
                raise TypeError('Expected rank <=2 dense array or matrix.')

            if arg1.ndim == 1:
                if dtype is not None:
                    arg1 = arg1.astype(dtype)
                sd.update((i, v) for i, v in enumerate(arg1.ravel()) if v != 0)
                self.dtype = arg1.dtype
            else:
                d = self._coo_container(arg1, dtype=dtype).todok()
                sd.update(d._dict.items())
                self.dtype = d.dtype
            self._shape = check_shape(arg1.shape, allow_1d=is_array)

    def update(self, val):
        # Prevent direct usage of update
        raise NotImplementedError(
            "Direct modification to dok_array element is not allowed."
        )

    def _update(self, data):
        """An update method for dict data defined for direct access to
        `dok_array` data. Main purpose is to be used for efficient conversion
        from other _spbase classes. Has no checking if `data` is valid.
        Be careful to use int for 1d-array and tuple for 2d.
        """
        return self._dict.update(data)

    def _getnnz(self, axis=None):
        if axis is not None:
            raise NotImplementedError(
                "_getnnz over an axis is not implemented for DOK format."
            )
        return len(self._dict)

    def count_nonzero(self):
        return sum(x != 0 for x in self._dict.values())

    _getnnz.__doc__ = _spbase._getnnz.__doc__
    count_nonzero.__doc__ = _spbase.count_nonzero.__doc__

    def __len__(self):
        return len(self._dict)

    def __contains__(self, key):
        return key in self._dict

    def setdefault(self, key, default=None, /):
        return self._dict.setdefault(key, default)

    def __delitem__(self, key, /):
        del self._dict[key]

    def clear(self):
        return self._dict.clear()

    def popitem(self):
        return self._dict.popitem()

    def items(self):
        return self._dict.items()

    def keys(self):
        return self._dict.keys()

    def values(self):
        return self._dict.values()

    def get(self, key, default=0.0):
        """This provides dict.get method functionality with type checking"""
        if key in self._dict:
            return self._dict.get(key, default)
        if isintlike(key) and self.ndim == 1:
            key = (key,)
        if self.ndim != len(key):
            raise IndexError(f'Index {key} length needs to match self.shape')
        try:
            for i in key:
                assert isintlike(i)
        except (AssertionError, TypeError, ValueError) as e:
            raise IndexError('Index must be or consist of integers.') from e
        key = tuple(i + M if i < 0 else i for i, M in zip(key, self.shape))
        if any(i < 0 or i >= M for i, M in zip(key, self.shape)):
            raise IndexError('Index out of bounds.')
        return self._dict.get(key, default)

    # 1D get methods
    def _get_int(self, idx):
        return self._dict.get(idx, self.dtype.type(0))

    def _get_slice(self, idx):
        i_range = range(*idx.indices(self.shape[0]))
        return self._get_array(list(i_range))

    def _get_array(self, idx):
        idx = np.asarray(idx)
        if len(idx.shape) == 0:
            val = self._dict.get(int(idx), self.dtype.type(0))
            return np.array(val, stype=self.dtype)
        new_dok = self._dok_container(idx.shape, dtype=self.dtype)
        dok_vals = [self._dict.get(i, 0) for i in idx.ravel()]
        if dok_vals:
            if len(idx.shape) == 1:
                for i, v in enumerate(dok_vals):
                    if v:
                        new_dok._dict[i] = v
            else:
                new_idx = np.unravel_index(np.arange(len(dok_vals)), idx.shape)
                new_idx = new_idx[0] if len(new_idx) == 1 else zip(*new_idx)
                # zip could use 'strict=True' With Python 3.10+
                for i, v in zip(new_idx, dok_vals):  # , strict=True):
                    if v:
                        new_dok._dict[i] = v
        return new_dok

    # 2D get methods
    def _get_intXint(self, row, col):
        return self._dict.get((row, col), self.dtype.type(0))

    def _get_intXslice(self, row, col):
        return self._get_sliceXslice(slice(row, row + 1), col)

    def _get_sliceXint(self, row, col):
        return self._get_sliceXslice(row, slice(col, col + 1))

    def _get_sliceXslice(self, row, col):
        row_start, row_stop, row_step = row.indices(self.shape[0])
        col_start, col_stop, col_step = col.indices(self.shape[1])
        row_range = range(row_start, row_stop, row_step)
        col_range = range(col_start, col_stop, col_step)
        shape = (len(row_range), len(col_range))
        # Switch paths only when advantageous
        # (count the iterations in the loops, adjust for complexity)
        if len(self) >= 2 * shape[0] * shape[1]:
            # O(nr*nc) path: loop over <row x col>
            return self._get_columnXarray(row_range, col_range)
        # O(nnz) path: loop over entries of self
        newdok = self._dok_container(shape, dtype=self.dtype)
        for key in self.keys():
            i, ri = divmod(int(key[0]) - row_start, row_step)
            if ri != 0 or i < 0 or i >= shape[0]:
                continue
            j, rj = divmod(int(key[1]) - col_start, col_step)
            if rj != 0 or j < 0 or j >= shape[1]:
                continue
            newdok._dict[i, j] = self._dict[key]
        return newdok

    def _get_intXarray(self, row, col):
        col = col.squeeze()
        return self._get_columnXarray([row], col)

    def _get_arrayXint(self, row, col):
        row = row.squeeze()
        return self._get_columnXarray(row, [col])

    def _get_sliceXarray(self, row, col):
        row = list(range(*row.indices(self.shape[0])))
        return self._get_columnXarray(row, col)

    def _get_arrayXslice(self, row, col):
        col = list(range(*col.indices(self.shape[1])))
        return self._get_columnXarray(row, col)

    def _get_columnXarray(self, row, col):
        # outer indexing
        newdok = self._dok_container((len(row), len(col)), dtype=self.dtype)

        for i, r in enumerate(row):
            for j, c in enumerate(col):
                v = self._dict.get((r, c), 0)
                if v:
                    newdok._dict[i, j] = v
        return newdok

    def _get_arrayXarray(self, row, col):
        # inner indexing
        i, j = map(np.atleast_2d, np.broadcast_arrays(row, col))
        newdok = self._dok_container(i.shape, dtype=self.dtype)

        for key in itertools.product(range(i.shape[0]), range(i.shape[1])):
            v = self._dict.get((i[key], j[key]), 0)
            if v:
                newdok._dict[key] = v
        return newdok

    # 1D set methods
    def _set_int(self, idx, x):
        if x:
            self._dict[idx] = x
        elif idx in self._dict:
            del self._dict[idx]

    def _set_slice(self, idx, x):
        i_range = range(*idx.indices(self.shape[0]))
        x = x.ravel()
        for i, v in zip(i_range, x):
            if v:
                self._dict[i] = v
            elif i in self._dict:
                del self._dict[i]

    def _set_array(self, idx, x):
        idx_set = idx.ravel()
        x_set = x.ravel()
        if len(idx_set) != len(x_set):
            if len(x_set) == 1:
                x_set = np.array([x_set[0]] * len(idx_set), dtype=self.dtype)
            else:
                raise ValueError("Need len(index)==len(data) or len(data)==1")
        for i, v in zip(idx_set, x_set):
            if v:
                self._dict[i] = v
            elif i in self._dict:
                del self._dict[i]

    # 2D set methods
    def _set_intXint(self, row, col, x):
        key = (row, col)
        if x:
            self._dict[key] = x
        elif key in self._dict:
            del self._dict[key]

    def _set_arrayXarray(self, row, col, x):
        row = list(map(int, row.ravel()))
        col = list(map(int, col.ravel()))
        x = x.ravel()
        self._dict.update(zip(zip(row, col), x))

        for i in np.nonzero(x == 0)[0]:
            key = (row[i], col[i])
            if self._dict[key] == 0:
                # may have been superseded by later update
                del self._dict[key]

    def __add__(self, other):
        if isscalarlike(other):
            res_dtype = upcast_scalar(self.dtype, other)
            new = self._dok_container(self.shape, dtype=res_dtype)
            # Add this scalar to every element.
            if self.ndim == 1:
                keys = range(self.shape[0])
            else:
                keys = itertools.product(*map(range, self.shape))
            for key in keys:
                aij = self._dict.get(key, 0) + other
                if aij:
                    new[key] = aij
        elif issparse(other):
            if other.shape != self.shape:
                raise ValueError("Matrix dimensions are not equal.")
            res_dtype = upcast(self.dtype, other.dtype)
            new = self._dok_container(self.shape, dtype=res_dtype)
            new._dict.update(self._dict)
            if other.format == "dok":
                with np.errstate(over='ignore'):
                    new._dict.update(
                        (k, new[k] + o_value) for k, o_value in other.items()
                    )
            else:
                o_coo = other.tocoo()
                if other.ndim == 1:
                    o_coo_data = zip(o_coo.indices[0], o_coo.data)
                else:
                    o_coo_data = zip(zip(*o_coo.indices), o_coo.data)
                with np.errstate(over='ignore'):
                    new._dict.update((k, new[k] + o_value) for k, o_value in o_coo_data)
        elif isdense(other):
            new = self.todense() + other
        else:
            return NotImplemented
        return new

    def __radd__(self, other):
        return self + other  # addition is commutative

    def __neg__(self):
        if self.dtype.kind == 'b':
            raise NotImplementedError(
                'Negating a sparse boolean matrix is not supported.'
            )
        new = self._dok_container(self.shape, dtype=self.dtype)
        new._dict.update((k, -v) for k, v in self.items())
        return new

    def _mul_scalar(self, other):
        res_dtype = upcast_scalar(self.dtype, other)
        # Multiply this scalar by every element.
        new = self._dok_container(self.shape, dtype=res_dtype)
        new._dict.update(((k, v * other) for k, v in self.items()))
        return new

    def _mul_vector(self, other):
        sd = self._dict
        res_dtype = upcast(self.dtype, other.dtype)
        # vector * vector
        if self.ndim == 1:
            if issparse(other):
                if other.format == "dok":
                    shared_keys = sd.keys() & other._dict.keys()
                else:
                    o_coo = other.tocoo()
                    shared_keys = sd.keys() & o_coo.indices[0]
                return np.array([sd[k] * other._dict[k] for k in shared_keys]).sum(
                    dtype=res_dtype
                )
            elif isdense(other):
                return np.array([other[k] * v for k, v in sd.items()]).sum(
                    dtype=res_dtype
                )

        # matrix * vector
        result = np.zeros(self.shape[0], dtype=res_dtype)
        for (i, j), v in sd.items():
            result[i] += v * other[j]
        return result

    def _mul_multivector(self, other):
        result_dtype = upcast(self.dtype, other.dtype)
        # vector * multivector
        if self.ndim == 1:
            if other.ndim == 1:
                result_shape = ()
                result = np.zeros(result_shape, dtype=result_dtype)
                for j, v in self._dict.items():
                    result += v * other[j]
            else:
                result_shape = (other.shape[1],)
                result = np.zeros(result_shape, dtype=result_dtype)
                for j, v in self._dict.items():
                    result[:] += v * other[j, :]
            return result

        # matrix * multivector
        if other.ndim == 1:
            result_shape = (self.shape[0],)
            result = np.zeros(result_shape, dtype=result_dtype)
            for (i, j), v in self.items():
                result[i] += v * other[j]
        else:
            result_shape = (self.shape[0], other.shape[1])
            result = np.zeros(result_shape, dtype=result_dtype)
            for (i, j), v in self.items():
                result[i, :] += v * other[j, :]
        return result

    def __imul__(self, other):
        if isscalarlike(other):
            self._dict.update((k, v * other) for k, v in self.items())
            return self
        return NotImplemented

    def __truediv__(self, other):
        if isscalarlike(other):
            res_dtype = upcast_scalar(self.dtype, other)
            new = self._dok_container(self.shape, dtype=res_dtype)
            new._dict.update(((k, v / other) for k, v in self._dict.items()))
            return new
        return self.tocsr() / other

    def __itruediv__(self, other):
        if isscalarlike(other):
            self._dict.update((k, v / other) for k, v in self.items())
            return self
        return NotImplemented

    def __reduce__(self):
        # this approach is necessary because __setstate__ is called after
        # __setitem__ upon unpickling and since __init__ is not called there
        # is no shape attribute hence it is not possible to unpickle it.
        return dict.__reduce__(self)

    def transpose(self, axes=None, copy=False):
        if self.ndim == 1:
            return self.copy()

        if axes is not None and axes != (1, 0):
            raise ValueError(
                "Sparse matrices do not support "
                "an 'axes' parameter because swapping "
                "dimensions is the only logical permutation."
            )

        M, N = self.shape
        new = self._dok_container((N, M), dtype=self.dtype, copy=copy)
        new._dict.update((((right, left), val) for (left, right), val in self.items()))
        return new

    transpose.__doc__ = _spbase.transpose.__doc__

    def conjtransp(self):
        """Return the conjugate transpose."""
        if self.ndim == 1:
            raise NotImplementedError

        M, N = self.shape
        new = self._dok_container((N, M), dtype=self.dtype)
        new._dict.update(
            (((right, left), np.conj(val)) for (left, right), val in self.items())
        )
        return new

    def copy(self):
        new = self._dok_container(self.shape, dtype=self.dtype)
        new._dict.update(self._dict)
        return new

    copy.__doc__ = _spbase.copy.__doc__

    def tocoo(self, copy=False):
        nnz = self.nnz
        sd = self._dict
        if nnz == 0:
            return self._coo_container(self.shape, dtype=self.dtype)

        idx_dtype = self._get_index_dtype(maxval=max(self.shape))
        data = np.fromiter(sd.values(), dtype=self.dtype, count=nnz)
        if self.ndim == 1:
            indices = (np.fromiter(sd, dtype=idx_dtype, count=nnz),)
        else:
            row = np.fromiter((i for i, _ in sd), dtype=idx_dtype, count=nnz)
            col = np.fromiter((j for _, j in sd), dtype=idx_dtype, count=nnz)
            indices = (row, col)
        A = self._coo_container((data, indices), shape=self.shape, dtype=self.dtype)
        A.has_canonical_format = True
        return A

    tocoo.__doc__ = _spbase.tocoo.__doc__

    def todok(self, copy=False):
        if copy:
            return self.copy()
        return self

    todok.__doc__ = _spbase.todok.__doc__

    def tocsc(self, copy=False):
        if self.ndim == 1:
            raise NotImplementedError("tocsr() not valid for 1d sparse array")
        return self.tocoo(copy=False).tocsc(copy=copy)

    tocsc.__doc__ = _spbase.tocsc.__doc__

    def eliminate_zeros(self):
        self._dict = {k: v for k, v in self._dict.items() if v != 0}

    def resize(self, *shape):
        is_array = isinstance(self, sparray)
        shape = check_shape(shape, allow_1d=is_array)
        if len(shape) != len(self.shape):
            # TODO implement resize across dimensions
            raise NotImplementedError

        if self.ndim == 1:
            newM = shape[0]
            for i in list(self._dict):
                if i >= newM:
                    del self._dict[i]
            self._shape = shape
            return

        newM, newN = shape
        M, N = self.shape
        if newM < M or newN < N:
            # Remove all elements outside new dimensions
            for i, j in list(self.keys()):
                if i >= newM or j >= newN:
                    del self._dict[i, j]
        self._shape = shape

    resize.__doc__ = _spbase.resize.__doc__

    def astype(self, dtype, casting='unsafe', copy=True):
        dtype = np.dtype(dtype)
        if self.dtype != dtype:
            result = self._dok_container(self.shape, dtype=dtype)
            data = np.array(list(self._dict.values()), dtype=dtype)
            result._update(zip(self._dict, data))
            return result
        elif copy:
            return self.copy()
        return self


def isspmatrix_dok(x):
    """Is `x` of dok_array type?

    Parameters
    ----------
    x
        object to check for being a dok matrix

    Returns
    -------
    bool
        True if `x` is a dok matrix, False otherwise

    Examples
    --------
    >>> from scipy.sparse import dok_array, dok_matrix, coo_matrix, isspmatrix_dok
    >>> isspmatrix_dok(dok_matrix([[5]]))
    True
    >>> isspmatrix_dok(dok_array([[5]]))
    False
    >>> isspmatrix_dok(coo_matrix([[5]]))
    False
    """
    return isinstance(x, dok_matrix)


# This namespace class separates array from matrix with isinstance
class dok_array(_dok_base, sparray):
    """
    Dictionary Of Keys based sparse array.

    This is an efficient structure for constructing sparse
    arrays incrementally.

    This can be instantiated in several ways:
        dok_array(D)
            where D is a 2-D ndarray

        dok_array(S)
            with another sparse array or matrix S (equivalent to S.todok())

        dok_array((M,N), [dtype])
            create the array with initial shape (M,N)
            dtype is optional, defaulting to dtype='d'

    Attributes
    ----------
    dtype : dtype
        Data type of the array
    shape : 2-tuple
        Shape of the array
    ndim : int
        Number of dimensions (this is always 2)
    nnz
        Number of nonzero elements
    size
    T

    Notes
    -----

    Sparse arrays can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    - Allows for efficient O(1) access of individual elements.
    - Duplicates are not allowed.
    - Can be efficiently converted to a coo_array once constructed.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import dok_array
    >>> S = dok_array((5, 5), dtype=np.float32)
    >>> for i in range(5):
    ...     for j in range(5):
    ...         S[i, j] = i + j    # Update element

    """


class dok_matrix(spmatrix, _dok_base, dict):
    """
    Dictionary Of Keys based sparse matrix.

    This is an efficient structure for constructing sparse
    matrices incrementally.

    This can be instantiated in several ways:
        dok_matrix(D)
            where D is a 2-D ndarray

        dok_matrix(S)
            with another sparse array or matrix S (equivalent to S.todok())

        dok_matrix((M,N), [dtype])
            create the matrix with initial shape (M,N)
            dtype is optional, defaulting to dtype='d'

    Attributes
    ----------
    dtype : dtype
        Data type of the matrix
    shape : 2-tuple
        Shape of the matrix
    ndim : int
        Number of dimensions (this is always 2)
    nnz
        Number of nonzero elements
    size
    T

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    - Allows for efficient O(1) access of individual elements.
    - Duplicates are not allowed.
    - Can be efficiently converted to a coo_matrix once constructed.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import dok_matrix
    >>> S = dok_matrix((5, 5), dtype=np.float32)
    >>> for i in range(5):
    ...     for j in range(5):
    ...         S[i, j] = i + j    # Update element

    """

    def set_shape(self, shape):
        new_matrix = self.reshape(shape, copy=False).asformat(self.format)
        self.__dict__ = new_matrix.__dict__

    def get_shape(self):
        """Get shape of a sparse matrix."""
        return self._shape

    shape = property(fget=get_shape, fset=set_shape)
