"""Dense 'dps' format for sparse array class. Dense Pretending Sparse"""
import numpy as np

from ._data import _data_matrix
from ._base import _sparray
from ._sputils import isshape, getdtype, check_reshape_kwargs, check_shape

__all__ = ["dps_array", "issparray_dps"]


MAXPRINT = 50


class dps_array(_data_matrix):
    format = "dps"
    _is_array = True

    def _with_data(self, data, copy=True):
        if copy:
            return self.__class__(data.copy(), shape=self.shape, dtype=data.dtype)
        else:
            return self.__class__(data, shape=self.shape, dtype=data.dtype)

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        _data_matrix.__init__(self)

        if isinstance(arg1, tuple):
            if isshape(arg1):
                M, N = arg1
                self._shape = check_shape((M, N))
                data_dtype = getdtype(dtype, default=float)
                self.data = np.zeros((M, N), dtype=data_dtype)
                self.has_canonical_format = True
            else:
                raise TypeError("invalid input format")
        else:
            M = np.atleast_2d(np.asarray(arg1))
            self._shape = check_shape(M.shape)
            self.data = M

        if dtype is not None:
            self.data = self.data.astype(dtype, copy=False)

        self._check()

    def reshape(self, *args, **kwargs):
        shape = check_shape(args, self.shape)
        order, copy = check_reshape_kwargs(kwargs)

        # Return early if reshape is not required
        if shape == self.shape:
            if copy:
                return self.copy()
            else:
                return self

        new_rows, new_cols = self.shape

        # Handle copy here rather than passing on to the constructor so that no
        # copy will be made of new_row and new_col regardless
        if copy:
            new_data = self.data.copy()
        else:
            new_data = self.data

        return self.__class__(new_data, shape=shape, copy=False)

    def getnnz(self, axis=None):
        return self.data.size

    getnnz.__doc__ = _sparray.getnnz.__doc__

    def _check(self):
        return True

    def transpose(self, axes=None, copy=False):
        if axes is not None:
            raise ValueError(
                "Sparse matrices do not support an 'axes' parameter because"
                "swapping dimensions is the only reasonable permutation."
            )

        M, N = self.shape
        return self.__class__(self.data.transpose(), copy=copy)

    transpose.__doc__ = _sparray.transpose.__doc__

    def resize(self, *shape):
        shape = check_shape(shape)
        new_M, new_N = shape
        M, N = self.shape

        if new_M < M or new_N < N:
            mask = np.logical_and(self.row < new_M, self.col < new_N)
            if not mask.all():
                self.data = self.data[mask]

        self._shape = shape

    resize.__doc__ = _sparray.resize.__doc__

    def toarray(self, order=None, out=None):
        return self.data.copy()

    def tocsc(self, copy=False):
        return self._csc_container(self.data)

    def tocsr(self, copy=False):
        return self._csr_container(self.data)

    def tocoo(self, copy=False):
        return self._coo_container(self.data)

    def todok(self, copy=False):
        return self._dok_container(self.data)

    def _add_dense(self, other):
        if other.shape != self.shape:
            raise ValueError(
                "Incompatible shapes ({} and {})".format(self.shape, other.shape)
            )
        return self._container(self.data, copy=False)

    def _mul_vector(self, other):
        return self.data @ other

    def _mul_multivector(self, other):
        return self.data @ other


def issparray_dps(x):
    return isinstance(x, dps_array)
