"""Compressed Sparse Row matrix format"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['csr_matrix', 'isspmatrix_csr']


from warnings import warn

import numpy as np
from scipy.lib.six import xrange

from .sparsetools import csr_tocsc, csr_tobsr, csr_count_blocks, \
        get_csr_submatrix, csr_sample_values
from .sputils import upcast, isintlike, IndexMixin, issequence

from .compressed import _cs_matrix


class csr_matrix(_cs_matrix, IndexMixin):
    """
    Compressed Sparse Row matrix

    This can be instantiated in several ways:
        csr_matrix(D)
            with a dense matrix or rank-2 ndarray D

        csr_matrix(S)
            with another sparse matrix S (equivalent to S.tocsr())

        csr_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        csr_matrix((data, ij), [shape=(M, N)])
            where ``data`` and ``ij`` satisfy the relationship
            ``a[ij[0, k], ij[1, k]] = data[k]``

        csr_matrix((data, indices, indptr), [shape=(M, N)])
            is the standard CSR representation where the column indices for
            row i are stored in ``indices[indptr[i]:indptr[i+1]]`` and their
            corresponding values are stored in ``data[indptr[i]:indptr[i+1]]``.
            If the shape parameter is not supplied, the matrix dimensions
            are inferred from the index arrays.

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
    data
        CSR format data array of the matrix
    indices
        CSR format index array of the matrix
    indptr
        CSR format index pointer array of the matrix
    has_sorted_indices
        Whether indices are sorted

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Advantages of the CSR format
      - efficient arithmetic operations CSR + CSR, CSR * CSR, etc.
      - efficient row slicing
      - fast matrix vector products

    Disadvantages of the CSR format
      - slow column slicing operations (consider CSC)
      - changes to the sparsity structure are expensive (consider LIL or DOK)

    Examples
    --------

    >>> from scipy.sparse import *
    >>> from scipy import *
    >>> csr_matrix( (3,4), dtype=int8 ).todense()
    matrix([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]], dtype=int8)

    >>> row = array([0,0,1,2,2,2])
    >>> col = array([0,2,2,0,1,2])
    >>> data = array([1,2,3,4,5,6])
    >>> csr_matrix( (data,(row,col)), shape=(3,3) ).todense()
    matrix([[1, 0, 2],
            [0, 0, 3],
            [4, 5, 6]])

    >>> indptr = array([0,2,3,6])
    >>> indices = array([0,2,2,0,1,2])
    >>> data = array([1,2,3,4,5,6])
    >>> csr_matrix( (data,indices,indptr), shape=(3,3) ).todense()
    matrix([[1, 0, 2],
            [0, 0, 3],
            [4, 5, 6]])

    """

    def transpose(self, copy=False):
        from .csc import csc_matrix
        M,N = self.shape
        return csc_matrix((self.data,self.indices,self.indptr), shape=(N,M), copy=copy)

    def tolil(self):
        from .lil import lil_matrix
        lil = lil_matrix(self.shape,dtype=self.dtype)

        self.sort_indices()  # lil_matrix needs sorted column indices

        ptr,ind,dat = self.indptr,self.indices,self.data
        rows, data = lil.rows, lil.data

        for n in xrange(self.shape[0]):
            start = ptr[n]
            end = ptr[n+1]
            rows[n] = ind[start:end].tolist()
            data[n] = dat[start:end].tolist()

        return lil

    def tocsr(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsc(self):
        indptr = np.empty(self.shape[1] + 1, dtype=np.intc)
        indices = np.empty(self.nnz, dtype=np.intc)
        data = np.empty(self.nnz, dtype=upcast(self.dtype))

        csr_tocsc(self.shape[0], self.shape[1],
                  self.indptr, self.indices, self.data,
                  indptr, indices, data)

        from .csc import csc_matrix
        A = csc_matrix((data, indices, indptr), shape=self.shape)
        A.has_sorted_indices = True
        return A

    def tobsr(self, blocksize=None, copy=True):
        from .bsr import bsr_matrix

        if blocksize is None:
            from .spfuncs import estimate_blocksize
            return self.tobsr(blocksize=estimate_blocksize(self))

        elif blocksize == (1,1):
            arg1 = (self.data.reshape(-1,1,1),self.indices,self.indptr)
            return bsr_matrix(arg1, shape=self.shape, copy=copy)

        else:
            R,C = blocksize
            M,N = self.shape

            if R < 1 or C < 1 or M % R != 0 or N % C != 0:
                raise ValueError('invalid blocksize %s' % blocksize)

            blks = csr_count_blocks(M,N,R,C,self.indptr,self.indices)

            indptr = np.empty(M//R + 1, dtype=np.intc)
            indices = np.empty(blks, dtype=np.intc)
            data = np.zeros((blks,R,C), dtype=self.dtype)

            csr_tobsr(M, N, R, C, self.indptr, self.indices, self.data,
                    indptr, indices, data.ravel())

            return bsr_matrix((data,indices,indptr), shape=self.shape)

    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[0],x[1])

    def __getitem__(self, key):
        def asindices(x):
            try:
                x = np.asarray(x, dtype=np.intc)
            except:
                raise IndexError('invalid index')
            else:
                return x

        def check_bounds(indices,N):
            max_indx = indices.max()
            if max_indx >= N:
                raise IndexError('index (%d) out of range' % max_indx)

            min_indx = indices.min()
            if min_indx < -N:
                raise IndexError('index (%d) out of range' % (N + min_indx))

            return (min_indx,max_indx)

        def extractor(indices,N):
            """Return a sparse matrix P so that P*self implements
            slicing of the form self[[1,2,3],:]
            """
            indices = asindices(indices)

            (min_indx,max_indx) = check_bounds(indices,N)

            if min_indx < 0:
                indices = indices.copy()
                indices[indices < 0] += N

            indptr = np.arange(len(indices) + 1, dtype=np.intc)
            data = np.ones(len(indices), dtype=self.dtype)
            shape = (len(indices),N)

            return csr_matrix((data,indices,indptr), shape=shape)

        row, col = self._unpack_index(key)

        # First attempt to use original row optimized methods
        # [1, ?]
        if isintlike(row):
            # [i, j]
            if isintlike(col):
                return self._get_single_element(row, col)
            # [i, 1:2]
            elif isinstance(col, slice):
                return self._get_row_slice(row, col)
            # [i, [1, 2]]
            elif issequence(col):
                P = extractor(col,self.shape[1]).T
                return self[row, :] * P
        elif isinstance(row, slice):
            # [1:2,??]
            if ((isintlike(col) and row.step in (1, None)) or
                    (isinstance(col, slice) and
                     col.step in (1, None) and
                     row.step in (1, None))):
                # col is int or slice with step 1, row is slice with step 1.
                return self._get_submatrix(row, col)
            elif issequence(col):
                P = extractor(col,self.shape[1]).T        # [1:2,[1,2]]
                # row is slice, col is sequence.
                return self[row,:]*P
        elif issequence(row):
            # [[1,2],??]
            if isintlike(col) or isinstance(col,slice):
                P = extractor(row, self.shape[0])     # [[1,2],j] or [[1,2],1:2]
                return (P*self)[:,col]
            elif issequence(col):                     # [[1,2],[1,2]]
                row = asindices(row)
                col = asindices(col)
                if row.ndim == 1:
                    if row.shape != col.shape:
                        raise IndexError('number of row and column indices differ')
                    check_bounds(row, self.shape[0])
                    check_bounds(col, self.shape[1])

                    num_samples = len(row)
                    val = np.empty(num_samples, dtype=self.dtype)
                    csr_sample_values(self.shape[0], self.shape[1],
                                      self.indptr, self.indices, self.data,
                                      num_samples, row, col, val)
                    return np.asmatrix(val)

        # If all else fails, try elementwise
        row, col = self._index_to_arrays(row, col)

        return self.__class__([[self._get_single_element(iii, jjj) for
                                iii, jjj in zip(ii, jj)] for ii, jj in
                               zip(row.tolist(), col.tolist())])

    def _get_single_element(self,row,col):
        """Returns the single element self[row, col]
        """
        M, N = self.shape
        if (row < 0):
            row += M
        if (col < 0):
            col += N
        if not (0 <= row < M) or not (0 <= col < N):
            raise IndexError("index out of bounds")

        # TODO make use of sorted indices (if present)

        start = self.indptr[row]
        end = self.indptr[row+1]
        indxs = np.where(col == self.indices[start:end])[0]

        num_matches = len(indxs)

        if num_matches == 0:
            # entry does not appear in the matrix
            return self.dtype.type(0)
        elif num_matches == 1:
            return self.data[start:end][indxs[0]]
        else:
            raise ValueError('nonzero entry (%d,%d) occurs more than once' % (row,col))

    def getrow(self, i):
        """Returns a copy of row i of the matrix, as a (1 x n)
        CSR matrix (row vector).
        """
        return self._get_submatrix(i, slice(None))

    def getcol(self, i):
        """Returns a copy of column i of the matrix, as a (m x 1)
        CSR matrix (column vector).
        """
        return self._get_submatrix(slice(None), i)

    def _get_row_slice(self, i, cslice):
        """Returns a copy of row self[i, cslice]
        """
        if i < 0:
            i += self.shape[0]

        if i < 0 or i >= self.shape[0]:
            raise IndexError('index (%d) out of range' % i)

        start, stop, stride = cslice.indices(self.shape[1])

        if stride == 1:
            # for stride == 1, _get_submatrix is ~30% faster than below
            row_slice = self._get_submatrix(i, cslice)

        else:
            # other strides need new code
            row_indices = self.indices[self.indptr[i]:self.indptr[i + 1]]
            row_data = self.data[self.indptr[i]:self.indptr[i + 1]]

            if stride > 0:
                ind = (row_indices >= start) & (row_indices < stop)
            elif stride < 0:
                ind = (row_indices <= start) & (row_indices > stop)

            if abs(stride) > 1:
                ind = ind & ((row_indices - start) % stride == 0)

            row_indices = (row_indices[ind] - start) // stride
            row_data = row_data[ind]
            row_indptr = np.array([0, len(row_indices)])

            if stride < 0:
                row_data = row_data[::-1]
                row_indices = abs(row_indices[::-1])

            shape = (1, int(np.ceil(float(stop - start) / stride)))

            row_slice = csr_matrix((row_data, row_indices, row_indptr),
                                   shape=shape)

        return row_slice

    def _get_submatrix(self, row_slice, col_slice):
        """Return a submatrix of this matrix (new matrix is created)."""

        M,N = self.shape

        def process_slice(sl, num):
            if isinstance(sl, slice):
                if sl.step not in (1, None):
                    raise ValueError('slicing with step != 1 not supported')
                i0, i1 = sl.start, sl.stop
                if i0 is None:
                    i0 = 0
                elif i0 < 0:
                    i0 = num + i0

                if i1 is None:
                    i1 = num
                elif i1 < 0:
                    i1 = num + i1
                return i0, i1

            elif isintlike(sl):
                if sl < 0:
                    sl += num
                return sl, sl + 1
            else:
                raise TypeError('expected slice or scalar')

        def check_bounds(i0, i1, num):
            if not (0 <= i0 <= num) or not (0 <= i1 <= num) or not (i0 <= i1):
                raise IndexError(
                      "index out of bounds: 0 <= %d <= %d, 0 <= %d <= %d,"
                       " %d <= %d" %
                      (i0, num, i1, num, i0, i1))

        i0, i1 = process_slice(row_slice, M)
        j0, j1 = process_slice(col_slice, N)
        check_bounds(i0, i1, M)
        check_bounds(j0, j1, N)

        indptr, indices, data = get_csr_submatrix(M, N,
                self.indptr, self.indices, self.data,
                int(i0), int(i1), int(j0), int(j1))

        shape = (i1 - i0, j1 - j0)

        return self.__class__((data,indices,indptr), shape=shape)


def isspmatrix_csr(x):
    return isinstance(x, csr_matrix)
