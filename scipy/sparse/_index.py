"""Indexing mixin for sparse array/matrix classes.
"""
import numpy as np
from ._sputils import isintlike
from ._base import sparray, issparse

INT_TYPES = (int, np.integer)


def _broadcast_arrays(a, b):
    """
    Same as np.broadcast_arrays(a, b) but old writeability rules.

    NumPy >= 1.17.0 transitions broadcast_arrays to return
    read-only arrays. Set writeability explicitly to avoid warnings.
    Retain the old writeability rules, as our Cython code assumes
    the old behavior.
    """
    x, y = np.broadcast_arrays(a, b)
    x.flags.writeable = a.flags.writeable
    y.flags.writeable = b.flags.writeable
    return x, y


class IndexMixin:
    """
    This class provides common dispatching and validation logic for indexing.
    """
    def __getitem__(self, key):
        index, new_shape = self._validate_indices(key)

        # 1D array
        if len(index) == 1:
            idx = index[0]
            if isinstance(idx, np.ndarray):
                if idx.shape == ():
                    idx = idx.item()
            if isinstance(idx, INT_TYPES):
                res = self._get_int(idx)
            elif isinstance(idx, slice):
                res = self._get_slice(idx)
            else:  # assume array idx
                res = self._get_array(idx)

        else:  # 2D array
            ixtypes = tuple(
                0 if isinstance(ix, INT_TYPES) else 1 if isinstance(ix, slice) else 2
                for ix in index
            )

            if ixtypes[0] == 0:
                if ixtypes[1] == 0:
                    res = self._get_intXint(*index)
                elif ixtypes[1] == 1:
                    res = self._get_intXslice(*index)
                elif ixtypes[1] == 2:
                    res = self._get_intXarray(*index)
            elif ixtypes[0] == 1:
                if ixtypes[1] == 0:
                    res = self._get_sliceXint(*index)
                if ixtypes[1] == 1:
                    res = self._get_sliceXslice(*index)
                if ixtypes[1] == 2:
                    res = self._get_sliceXarray(*index)
            elif ixtypes[0] == 2:
                row, col = index
                if ixtypes[1] == 0:
                    res = self._get_arrayXint(*index)
                elif ixtypes[1] == 1:
                    if row.ndim == 2:
                        assert False, "You should get stopped in validate_indices before this"
                        raise IndexError('index results in >2 dimensions')
                    res = self._get_arrayXslice(*index)
                else:
                    # arrayXarray preprocess
                    if row.ndim == 2 and row.shape[1] == 1\
                        and (col.ndim == 1 or col.shape[0] == 1):
                        # outer indexing
                        res = self._get_columnXarray(row[:, 0], col.ravel())
                    else:
                        # inner indexing
                        row, col = _broadcast_arrays(*(np.array(ix) for ix in index))
                        if row.shape != col.shape:
                            raise IndexError('number of row and column indices differ')
                        if row.size == 0:
                            row_shape = np.atleast_2d(row).shape
                            res = self.__class__(row_shape, dtype=self.dtype)
                        else:
                            res = self._get_arrayXarray(row, col)

        if not isinstance(self, sparray):
            return res
        if self.format == 'lil' and len(new_shape) != 2:
            return res.tocoo().reshape(new_shape)
        if res.shape == () and new_shape != ():
            return self.__class__([res], shape=new_shape, dtype=self.dtype)
        return res.reshape(new_shape)

    def __setitem__(self, key, x):
        index, new_shape = self._validate_indices(key)

        # 1D array
        if len(index) == 1:
            idx = index[0]

            if issparse(x):
                x = x.toarray()
            x = np.asarray(x, dtype=self.dtype)

            if isinstance(idx, INT_TYPES):
                if x.size != 1:
                    raise ValueError('Trying to assign a sequence to an item')
                self._set_int(idx, x.flat[0])
                return

            if isinstance(idx, slice):
                # Note: Python `range` does not blow up memory even if shape is large
                idx_range = range(*idx.indices(self.shape[0]))
                N = len(idx_range)
                if N == 1 and x.size == 1:
                    self._set_int(idx_range[0], x.flat[0])

                # broadcast scalar to full 1d
                if x.squeeze().shape != (N,):
                    x = np.broadcast_to(x, (N,))
                if x.size != N:
                    raise ValueError(f'size mismatch: put {x.size} values in {N} spots')
                if x.size == 0:
                    return
                self._set_slice(idx, x)
                return

            # broadcast scalar to full 1d
            if x.squeeze().shape != idx.squeeze().shape:
                x = np.broadcast_to(x, idx.shape)
            if x.size == 0:
                return
            self._set_array(idx, x)
            return

        # 2D array

        row, col = index

        if isinstance(row, INT_TYPES) and isinstance(col, INT_TYPES):
            x = np.asarray(x, dtype=self.dtype)
            if x.size != 1:
                raise ValueError('Trying to assign a sequence to an item')
            self._set_intXint(row, col, x.flat[0])
            return

        if isinstance(row, slice):
            row = np.arange(*row.indices(self.shape[0]))[:, None]
        else:
            row = np.atleast_1d(row)

        if isinstance(col, slice):
            col = np.arange(*col.indices(self.shape[1]))[None, :]
            if row.ndim == 1:
                row = row[:, None]
        else:
            col = np.atleast_1d(col)

        i, j = _broadcast_arrays(row, col)
        if i.shape != j.shape:
            raise IndexError('number of row and column indices differ')

        if issparse(x):
            if i.ndim == 1:
                # Inner indexing, so treat them like row vectors.
                i = i[None]
                j = j[None]
            broadcast_row = x.shape[0] == 1 and i.shape[0] != 1
            broadcast_col = x.shape[1] == 1 and i.shape[1] != 1
            if not ((broadcast_row or x.shape[0] == i.shape[0]) and
                    (broadcast_col or x.shape[1] == i.shape[1])):
                raise ValueError('shape mismatch in assignment')
            if x.shape[0] == 0 or x.shape[1] == 0:
                return
            x = x.tocoo(copy=True)
            x.sum_duplicates()
            self._set_arrayXarray_sparse(i, j, x)
        else:
            # Make x and i into the same shape
            x = np.asarray(x, dtype=self.dtype)
            if x.squeeze().shape != i.squeeze().shape:
                x = np.broadcast_to(x, i.shape)
            if x.size == 0:
                return
            x = x.reshape(i.shape)
            self._set_arrayXarray(i, j, x)

    def _validate_indices(self, key):
        """Returns index tuple of intended result"""
        # single ellipsis
        if key is Ellipsis:
            return (slice(None),) * self.ndim, self.shape

        if not isinstance(key, tuple):
            key = [key]

        ellps_pos = None
        idx_shape = []
        index = []
        index_ndim = 0
        array_indices = []
        for i, idx in enumerate(key):
            if idx is Ellipsis:
                if ellps_pos is not None:
                    raise IndexError('an index can only have a single ellipsis')
                ellps_pos = i
            elif idx is None:
                idx_shape.append(1)
            elif isinstance(idx, slice):
                index.append(idx)
                len_slice = len(range(*idx.indices(self._shape[index_ndim])))
                idx_shape.append(len_slice)
                index_ndim += 1
            elif isintlike(idx):
                N = self._shape[index_ndim]
                if not (-N <= idx < N):
                    raise IndexError(f'index ({idx}) out of range')
                idx = int(idx + N if idx < 0 else idx)
                index.append(idx)
                index_ndim += 1
            elif (ix := self._compatible_boolean_index(idx)) is not None:
                tmp_ndim = index_ndim + ix.ndim
                mid_shape = self._shape[index_ndim:tmp_ndim]
                if ix.shape != mid_shape:
                    raise IndexError(
                        f"bool index {i} has shape {mid_shape} instead of {ix.shape}"
                    )
                index.extend(ix.nonzero())
                array_indices.extend(range(index_ndim, tmp_ndim))
                index_ndim = tmp_ndim
            elif issparse(idx):
                # TODO: make sparse matrix indexing work for sparray
                raise IndexError(
                    'Indexing with sparse matrices is not supported '
                    'except boolean indexing where matrix and index '
                    'are equal shapes.')
            else:  # dense array
                N = self._shape[index_ndim]
                idx = np.array(idx)
                idx = self._asindices(idx, N)
                index.append(idx)
                array_indices.append(index_ndim)
                index_ndim += 1
        if index_ndim > self.ndim:
            raise IndexError(
                f'invalid index ndim. Array is {self.ndim}D. Index needs {index_ndim}D'
            )
                # raise IndexError('invalid ndim for array index')
        if len(array_indices) > 1:
            idx_arrays = _broadcast_arrays(*(index[i] for i in array_indices))
            if any(idx_arrays[0].shape != ix.shape for ix in idx_arrays[1:]):
                raise IndexError('array indices after broadcast differ in shape')

        # add slice(None) (which is colon) to fill out full index
        nslice = self.ndim - index_ndim
        if nslice > 0:
            if ellps_pos is None:
                ellps_pos = index_ndim
            index = index[:ellps_pos] + [slice(None)] * nslice + index[ellps_pos:]
            mid_shape = list(self.shape[ellps_pos:ellps_pos + nslice])
            idx_shape = idx_shape[:ellps_pos] + mid_shape + idx_shape[ellps_pos:]

        if array_indices:
            idx_shape = list(index[array_indices[0]].shape) + idx_shape
        if len(idx_shape) > 2:
            ndim = len(idx_shape)
            raise IndexError(f'Only 1D or 2D arrays allowed. Index makes {ndim}D')
        return tuple(index), tuple(idx_shape)

    def _compatible_boolean_index(self, idx):
        """Check for boolean array or array-like. peek before asarray for array-like"""
        # assume already an array if attr ndim exists: skip to bottom
        if not hasattr(idx, 'ndim'):
            # is first element boolean?
            try:
                ix = next(iter(idx), None)
                for _ in range(self.ndim):
                    if isinstance(ix, bool):
                        break
                    ix = next(iter(ix), None)
                else:
                    return None
            except TypeError:
                return None
            # since first is boolean, construct array and check all elements
            idx = np.asanyarray(idx)

        if idx.dtype.kind == 'b':
            return idx

    def _asindices(self, idx, length):
        """Convert `idx` to a valid index for an axis with a given length.

        Subclasses that need special validation can override this method.
        """
        try:
            x = np.asarray(idx)
        except (ValueError, TypeError, MemoryError) as e:
            raise IndexError('invalid index') from e

        if x.ndim not in (1, 2):
            raise IndexError('Index dimension must be 1 or 2')

        if x.size == 0:
            return x

        # Check bounds
        max_indx = x.max()
        if max_indx >= length:
            raise IndexError('index (%d) out of range' % max_indx)

        min_indx = x.min()
        if min_indx < 0:
            if min_indx < -length:
                raise IndexError('index (%d) out of range' % min_indx)
            if x is idx or not x.flags.owndata:
                x = x.copy()
            x[x < 0] += length
        return x

    def _get_int(self, idx):
        raise NotImplementedError()

    def _get_slice(self, idx):
        raise NotImplementedError()

    def _get_array(self, idx):
        raise NotImplementedError()

    def _getrow(self, i):
        """Return a copy of row i of the matrix, as a (1 x n) row vector.
        """
        M, N = self.shape
        i = int(i)
        if i < -M or i >= M:
            raise IndexError('index (%d) out of range' % i)
        if i < 0:
            i += M
        return self._get_intXslice(i, slice(None))

    def _getcol(self, i):
        """Return a copy of column i of the matrix, as a (m x 1) column vector.
        """
        M, N = self.shape
        i = int(i)
        if i < -N or i >= N:
            raise IndexError('index (%d) out of range' % i)
        if i < 0:
            i += N
        return self._get_sliceXint(slice(None), i)

    def _get_intXint(self, row, col):
        raise NotImplementedError()

    def _get_intXarray(self, row, col):
        raise NotImplementedError()

    def _get_intXslice(self, row, col):
        raise NotImplementedError()

    def _get_sliceXint(self, row, col):
        raise NotImplementedError()

    def _get_sliceXslice(self, row, col):
        raise NotImplementedError()

    def _get_sliceXarray(self, row, col):
        raise NotImplementedError()

    def _get_arrayXint(self, row, col):
        raise NotImplementedError()

    def _get_arrayXslice(self, row, col):
        raise NotImplementedError()

    def _get_columnXarray(self, row, col):
        raise NotImplementedError()

    def _get_arrayXarray(self, row, col):
        raise NotImplementedError()

    def _set_int(self, idx, x):
        raise NotImplementedError()

    def _set_slice(self, idx, x):
        raise NotImplementedError()

    def _set_array(self, idx, x):
        raise NotImplementedError()

    def _set_intXint(self, row, col, x):
        raise NotImplementedError()

    def _set_arrayXarray(self, row, col, x):
        raise NotImplementedError()

    def _set_arrayXarray_sparse(self, row, col, x):
        # Fall back to densifying x
        x = np.asarray(x.toarray(), dtype=self.dtype)
        x, _ = _broadcast_arrays(x, row)
        self._set_arrayXarray(row, col, x)
