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
    def _raise_on_1d_array_slice(self):
        """We do not currently support 1D sparse arrays.

        This function is called each time that a 1D array would
        result, raising an error instead.

        Once 1D sparse arrays are implemented, it should be removed.
        """
        if isinstance(self, sparray):
            raise NotImplementedError(
                'We have not yet implemented 1D sparse slices; '
                'please index using explicit indices, e.g. `x[:, [0]]`'
            )

    def __getitem__(self, key):
        # handle 1d indexing
        if self.ndim == 1:
            idx = self._validate_indices(key)
            if isinstance(idx, tuple) and len(idx) == 1:
                idx = idx[0]
            if isinstance(idx, INT_TYPES):
                return self._get_int(idx)
            elif isinstance(idx, slice):
                return self._get_slice(idx)
            # assume array idx
            return self._get_array(idx)

        row, col = self._validate_indices(key)

        # Dispatch to specialized methods.
        if isinstance(row, INT_TYPES):
            if isinstance(col, INT_TYPES):
                return self._get_intXint(row, col)
            elif isinstance(col, slice):
                self._raise_on_1d_array_slice()
                return self._get_intXslice(row, col)
            elif col.ndim == 1:
                self._raise_on_1d_array_slice()
                return self._get_intXarray(row, col)
            elif col.ndim == 2:
                return self._get_intXarray(row, col)
            raise IndexError('index results in >2 dimensions')
        elif isinstance(row, slice):
            if isinstance(col, INT_TYPES):
                self._raise_on_1d_array_slice()
                return self._get_sliceXint(row, col)
            elif isinstance(col, slice):
                if row == slice(None) and row == col:
                    return self.copy()
                return self._get_sliceXslice(row, col)
            elif col.ndim == 1:
                return self._get_sliceXarray(row, col)
            raise IndexError('index results in >2 dimensions')
        elif row.ndim == 1:
            if isinstance(col, INT_TYPES):
                self._raise_on_1d_array_slice()
                return self._get_arrayXint(row, col)
            elif isinstance(col, slice):
                return self._get_arrayXslice(row, col)
        else:  # row.ndim == 2
            if isinstance(col, INT_TYPES):
                return self._get_arrayXint(row, col)
            elif isinstance(col, slice):
                raise IndexError('index results in >2 dimensions')
            elif row.shape[1] == 1 and (col.ndim == 1 or col.shape[0] == 1):
                # special case for outer indexing
                return self._get_columnXarray(row[:,0], col.ravel())

        # The only remaining case is inner (fancy) indexing
        row, col = _broadcast_arrays(row, col)
        if row.shape != col.shape:
            raise IndexError('number of row and column indices differ')
        if row.size == 0:
            return self.__class__(np.atleast_2d(row).shape, dtype=self.dtype)
        return self._get_arrayXarray(row, col)

    def __setitem__(self, key, x):
        # handle 1d indexing
        if self.ndim == 1:
            idx = self._validate_indices(key)
            if isinstance(idx, tuple) and len(idx) == 1:
                idx = idx[0]
            if isinstance(idx, INT_TYPES):
                x = np.asarray(x, dtype=self.dtype)
                if x.size != 1:
                    raise ValueError('Trying to assign a sequence to an item')
                self._set_int(idx, x.flat[0])
                return

            if isinstance(idx, slice):
                idx = np.arange(*idx.indices(self.shape[0]))
            else:
                idx = np.atleast_1d(idx)

            # broadcast scalar to full 1d
            if issparse(x):
                x = x.toarray()
            x = np.asarray(x, dtype=self.dtype)
            if x.squeeze().shape != idx.squeeze().shape:
                x = np.broadcast_to(x, idx.shape)
            if x.size == 0:
                return
            x = x.reshape(idx.shape)
            self._set_array(idx, x)
            return

        row, col = self._validate_indices(key)

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
            return (slice(None),) * self.ndim

        if not isinstance(key, tuple):
            key = [key]

        ellps_pos = None
        none_indices = []
        idx_shape = []
        index = []
        index_ndim = 0
        for i, idx in enumerate(key):
            if idx is Ellipsis:
                if ellps_pos is not None:
                    raise IndexError('an index can only have a single ellipsis')
                ellps_pos = i
            elif idx is None:
                none_indices.append(i)
            elif isinstance(idx, slice):
                index.append(idx)
                idx_shape.append(len(range(*idx.indices(self._shape[index_ndim]))))
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
                index_ndim = tmp_ndim
                idx_shape.append(len(index[-1]))
            elif issparse(idx):
                # TODO: make sparse matrix indexing work for sparray
                raise IndexError(
                    'Indexing with sparse matrices is not supported '
                    'except boolean indexing where matrix and index '
                    'are equal shapes.')
            else:  # dense array
                N = self._shape[index_ndim]
                idx = self._asindices(idx, N)
                index.append(idx)
                index_ndim += 1
                idx_shape.append(idx.size)
            if index_ndim > self.ndim:
                raise IndexError('invalid ndim for array index')

        # add slice(None) (which is colon) to fill out full index
        nslice = self.ndim - index_ndim
        if nslice > 0:
            if ellps_pos is None:
                ellps_pos = len(index)
            index = index[:ellps_pos] + [slice(None)] * nslice + index[ellps_pos:]
            idx_shape = idx_shape[:ellps_pos] + [1] * nslice + idx_shape[ellps_pos:]

        # slide `None` into index at none_positions from input key
        # (about 500 times faster than repeated index.insert(pos, None))
        if none_indices:
            new_index, old_dx = [], 0
            for nn, ix in enumerate(none_indices):
                dx = ix - nn
                new_index.extend(index[old_dx:dx])
                new_index.append(None)
                old_dx = dx
            new_index.extend(index[old_dx:])
            index = new_index
        return tuple(index) # , tuple(idx_shape)

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

    def _set_intXint(self, row, col, x):
        raise NotImplementedError()

    def _set_arrayXarray(self, row, col, x):
        raise NotImplementedError()

    def _set_arrayXarray_sparse(self, row, col, x):
        # Fall back to densifying x
        x = np.asarray(x.toarray(), dtype=self.dtype)
        x, _ = _broadcast_arrays(x, row)
        self._set_arrayXarray(row, col, x)
