"""Fast LInked List sparse matrix class
"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['fast_lil_matrix', 'isspmatrix_fast_lil']

import numpy as np

from scipy._lib.six import xrange
from .base import spmatrix, isspmatrix
from .sputils import (getdtype, isshape, isscalarlike, IndexMixin,
                      upcast_scalar, get_index_dtype, isintlike,
                      get_index_dtype)
from . import _fastlil
from .lil import isspmatrix_lil


class fast_lil_matrix(spmatrix, IndexMixin):
    """Row-based sparse matrix

    This is a structure for constructing sparse matrices incrementally.
    Note that inserting a single item can take linear time in the worst case;
    to construct a matrix efficiently, make sure the items are pre-sorted by
    index, per row.

    This can be instantiated in several ways:
        fast_lil_matrix(D)
            with a dense matrix or rank-2 ndarray D

        fast_lil_matrix(S)
            with another sparse matrix S (equivalent to S.tofastlil())

        fast_lil_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

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

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Advantages of the LIL format
        - supports flexible slicing
        - changes to the matrix sparsity structure are efficient

    Disadvantages of the LIL format
        - arithmetic operations LIL + LIL are slow (consider CSR or CSC)
        - slow column slicing (consider CSC)
        - slow matrix vector products (consider CSR or CSC)

    Intended Usage
        - LIL is a convenient format for constructing sparse matrices
        - once a matrix has been constructed, convert to CSR or
          CSC format for fast arithmetic and matrix vector operations
        - consider using the COO format when constructing large matrices

    Data Structure
        - An array (``self.rows``) of rows, each of which is a sorted
          list of column indices of non-zero elements.
        - The corresponding nonzero values are stored in similar
          fashion in ``self.data``.


    """
    format = 'fastlil'

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        spmatrix.__init__(self)
        self.dtype = getdtype(dtype, arg1, default=float)
        self._idx_dtype = None

        # First get the shape
        if isspmatrix(arg1):

            A = arg1

            if dtype is not None:
                A = A.astype(dtype)

            self.shape = A.shape
            self.dtype = A.dtype

            if isspmatrix_fast_lil(A):
                if copy:
                    A = A.copy()
                self._idx_dtype = A._idx_dtype
                self._matrix = A._matrix
            else:
                A = A.tocsr()
                self._idx_dtype = A.indices.dtype
                self._matrix = self._get_matrix()

                if A.data.dtype == np.bool:
                    self._matrix.fromcsr(A.indices, A.indptr, A.data.astype(np.uint8))
                else:
                    self._matrix.fromcsr(A.indices, A.indptr, A.data)

        elif isinstance(arg1, tuple):
            if isshape(arg1):
                if shape is not None:
                    raise ValueError('invalid use of shape parameter')
                M, N = arg1
                self.shape = (M, N)
                self._idx_dtype = np.dtype(get_index_dtype(maxval=max(*self.shape)))
                self._matrix = self._get_matrix()
            else:
                raise TypeError('unrecognized lil_matrix constructor usage')
        else:
            # assume A is dense
            try:
                A = np.asmatrix(arg1)
            except TypeError:
                raise TypeError('unsupported matrix type')
            else:
                from .csr import csr_matrix
                A = csr_matrix(A, dtype=dtype)

                self.shape = A.shape
                self.dtype = A.dtype
                self._idx_dtype = A.indices.dtype
                self._matrix = self._get_matrix()

                if A.data.dtype == np.bool:
                    self._matrix.fromcsr(A.indices, A.indptr, A.data.astype(np.uint8))
                else:
                    self._matrix.fromcsr(A.indices, A.indptr, A.data)

    def _get_matrix(self):

        rows, cols = self.shape

        return getattr(_fastlil,
                       'fast_lil_matrix_{}_{}'.format(
                           str(self._idx_dtype),
                           str(self.dtype)))(np.intp(rows),
                                             np.intp(cols))

    def _set(self, row, col, value):

        self._matrix.safe_set(self._idx_dtype(row),
                              self._idx_dtype(col),
                              self.dtype.type(value))

    def _from_lil(self, lil_matrix):

        for row_idx, (row_indices, row_data) in enumerate(zip(lil_matrix.rows,
                                                              lil_matrix.data)):
            for col_idx, value in zip(row_indices, row_data):
                self._set(row_idx, col_idx, value)

    def set_shape(self, shape):
        shape = tuple(shape)

        if len(shape) != 2:
            raise ValueError("Only two-dimensional sparse arrays "
                             "are supported.")
        try:
            shape = int(shape[0]), int(shape[1])  # floats, other weirdness
        except:
            raise TypeError('invalid shape')

        if not (shape[0] >= 0 and shape[1] >= 0):
            raise ValueError('invalid shape')

        if (self._shape != shape) and (self._shape is not None):
            try:
                self = self.reshape(shape)
            except NotImplementedError:
                raise NotImplementedError("Reshaping not implemented for %s." %
                                          self.__class__.__name__)
        self._shape = shape

    shape = property(fget=spmatrix.get_shape, fset=set_shape)

    def __iadd__(self, other):
        self[:, :] = self + other
        return self

    def __isub__(self, other):
        self[:, :] = self - other
        return self

    def __imul__(self, other):
        if isscalarlike(other):
            self[:, :] = self * other
            return self
        else:
            return NotImplemented

    def __itruediv__(self, other):
        if isscalarlike(other):
            self[:, :] = self / other
            return self
        else:
            return NotImplemented

    def getnnz(self, axis=None):
        """Get the count of explicitly-stored values (nonzeros)

        Parameters
        ----------
        axis : None, 0, or 1
            Select between the number of values across the whole matrix, in
            each column, or in each row.
        """

        if axis is None:
            return self._matrix.getnnz(None)
        if axis < 0:
            axis += 2

        if axis in (0, 1):
            return self._matrix.getnnz(axis)
        else:
            raise ValueError('axis out of bounds')

    def count_nonzero(self):
        return self._matrix.count_nonzero()

    getnnz.__doc__ = spmatrix.getnnz.__doc__
    count_nonzero.__doc__ = spmatrix.count_nonzero.__doc__

    nnz = property(fget=getnnz)

    def getrow(self, i):
        """Returns a copy of the 'i'th row.
        """

        return self[i, :]

    def __getitem__(self, index):
        """Return the element(s) index=(i, j), where j may be a slice.
        This always returns a copy for consistency, since slices into
        Python lists return copies.
        """

        # Scalar fast path first
        if isinstance(index, tuple) and len(index) == 2:
            i, j = index
            # Use isinstance checks for common index types; this is
            # ~25-50% faster than isscalarlike. Other types are
            # handled below.
            if ((isinstance(i, int) or isinstance(i, np.integer)) and
                    (isinstance(j, int) or isinstance(j, np.integer))):

                return self._matrix.safe_get(i, j)

        # Utilities found in IndexMixin
        i, j = self._unpack_index(index)

        i_intlike = False
        i_slice = False
        i_list = False

        j_intlike = False
        j_slice = False
        j_list = False

        # Proper check for other scalar index types
        if isintlike(i):
            i_intlike = True
        elif isinstance(i, slice):
            i_slice = True
        elif isinstance(i, (list, tuple, np.ndarray)):
            i_list = True
        else:
            raise ValueError

        if isintlike(j):
            j_intlike = True
        elif isinstance(j, slice):
            j_slice = True
        elif isinstance(j, (list, tuple, np.ndarray)):
            j_list = True
        else:
            raise ValueError

        # Fast path for integer indices
        if i_intlike and j_intlike:
            return self._matrix.safe_get(i, j)

        # Fast path for row indexing when we want all the cols
        # and the row index is an integer or 1d
        if j_slice and j == slice(None, None, None):
            if i_slice:
                i = np.arange(*i.indices(self.shape[0]),
                              dtype=self._matrix.idx_dtype())
            else:
                i = np.atleast_1d(i).astype(self._matrix.idx_dtype())

            if i.ndim == 1:
                new = fast_lil_matrix((i.shape[0],
                                       self.shape[1]),
                                      dtype=self.dtype)
                new._matrix = self._matrix.fancy_get_rows(i)

                return new

        # Fast path for col indexing when we want all the rows
        if i_slice and i == slice(None, None, None):
            if j_slice:
                j = np.arange(*j.indices(self.shape[1]),
                              dtype=self._matrix.idx_dtype())
            else:
                j = np.atleast_1d(j).astype(self._matrix.idx_dtype())

            if j.ndim == 1:
                new = fast_lil_matrix((self.shape[0],
                                       j.shape[0]),
                                      dtype=self.dtype)
                new._matrix = self._matrix.fancy_get_cols(j)

                return new

        # Full-blown indexing
        i, j = self._index_to_arrays(i, j)

        if i.size == 0:
            return fast_lil_matrix(i.shape, dtype=self.dtype)

        i, j = _prepare_index_for_memoryview(i, j, self._matrix.idx_dtype())
        new_rows, new_cols = i.shape

        new = fast_lil_matrix(i.shape, dtype=self.dtype)
        new._matrix = self._matrix.fancy_get(new_rows,
                                             new_cols,
                                             i,
                                             j)

        return new

    def __setitem__(self, index, x):
        # Scalar fast path first
        if isinstance(index, tuple) and len(index) == 2:
            i, j = index
            # Use isinstance checks for common index types; this is
            # ~25-50% faster than isscalarlike. Scalar index
            # assignment for other types is handled below together
            # with fancy indexing.
            if ((isinstance(i, int) or isinstance(i, np.integer)) and
                    (isinstance(j, int) or isinstance(j, np.integer))):
                x = self.dtype.type(x)
                if x.size > 1:
                    # Triggered if input was an ndarray
                    raise ValueError("Trying to assign a sequence to an item")

                self._matrix.safe_set(i, j, x)
                return

        # General indexing
        i, j = self._unpack_index(index)

        # shortcut for common case of full matrix assign:
        if (isspmatrix(x) and isinstance(i, slice) and i == slice(None) and
                isinstance(j, slice) and j == slice(None)
                and x.shape == self.shape):
            x = fast_lil_matrix(x, dtype=self.dtype)
            self._matrix = x._matrix
            return

        i, j = self._index_to_arrays(i, j)

        if isspmatrix(x):
            x = x.toarray()

        # Make x and i into the same shape
        x = np.asarray(x, dtype=self.dtype)
        x, _ = np.broadcast_arrays(x, i)

        if x.shape != i.shape:
            raise ValueError("shape mismatch in assignment")

        # Set values
        i, j, x = _prepare_index_for_memoryview(i, j,
                                                self._matrix.idx_dtype(), x=x)
        self._matrix.fancy_set(i, j, x)

    def _mul_scalar(self, other):
        if other == 0:
            # Multiply by zero: return the zero matrix
            new = fast_lil_matrix(self.shape, dtype=self.dtype)
        else:
            res_dtype = upcast_scalar(self.dtype, other)

            new = self.copy()
            new = new.astype(res_dtype)
            new._matrix.mul(new.dtype.type(other))

        return new

    def __truediv__(self, other):           # self / other
        if isscalarlike(other):
            new = self.copy()
            new._matrix.mul(1 / new.dtype.type(other))
            return new
        else:
            return self.tocsr() / other

    def copy(self):

        new = fast_lil_matrix(self.shape, dtype=self.dtype)
        new._matrix = self._matrix.copy()

        return new

    def reshape(self, shape):
        new = fast_lil_matrix(shape, dtype=self.dtype)

        new._matrix = self._matrix.reshape(*shape)

        return new

    def toarray(self, order=None, out=None):
        """See the docstring for `spmatrix.toarray`."""
        d = self._process_toarray_args(order, out)

        if d.dtype != np.bool:
            self._matrix.todense(d)
        else:
            for row_idx in range(self.shape[0]):
                row_indices, row_data = self._matrix.get_row(row_idx)

                for i in range(len(row_indices)):
                    d[row_idx, row_indices[i]] = row_data[i]

        return d

    def transpose(self):
        return self.tocsr().transpose().tocsr().tofastlil()

    def tofastlil(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self, copy=False):
        """ Return Compressed Sparse Row format arrays for this matrix.
        """

        indices, indptr, data = self._matrix.tocsr()

        if self.dtype == np.bool:
            data = data.astype(np.bool)

        from .csr import csr_matrix
        return csr_matrix((data, indices, indptr), shape=self.shape)

    def tocsc(self):
        """ Return Compressed Sparse Column format arrays for this matrix.
        """
        return self.tocsr().tocsc()

    def __getstate__(self):
        # Pickle through converting to and from CSR

        state = self.__dict__

        state['_matrix'] = self._matrix.tocsr()

        return state

    def __setstate__(self, state):

        self.__dict__.update(state)
        self._matrix = self._get_matrix()
        self._matrix.fromcsr(*state['_matrix'])


def _prepare_index_for_memoryview(i, j, idx_dtype, x=None):
    """
    Convert index and data arrays to form suitable for passing to the
    Cython fancy getset routines.

    The conversions are necessary since to (i) ensure the integer
    index arrays are in one of the accepted types, and (ii) to ensure
    the arrays are writable so that Cython memoryview support doesn't
    choke on them.

    Parameters
    ----------
    i, j
        Index arrays
    x : optional
        Data arrays

    Returns
    -------
    i, j, x
        Re-formatted arrays (x is omitted, if input was None)

    """

    if not i.flags.writeable or i.dtype is not idx_dtype:
        i = i.astype(idx_dtype)
    if not j.flags.writeable or j.dtype is not idx_dtype:
        j = j.astype(idx_dtype)

    if x is not None:

        if x.dtype == np.bool:
            x = x.astype(np.uint8)

        if not x.flags.writeable:
            x = x.copy()

        return i, j, x
    else:
        return i, j


def isspmatrix_fast_lil(x):
    return isinstance(x, fast_lil_matrix)
