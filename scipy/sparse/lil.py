"""LInked List sparse matrix class
"""

__docformat__ = "restructuredtext en"

__all__ = ['lil_matrix','isspmatrix_lil']

from bisect import bisect_left

import numpy as np

from base import spmatrix, isspmatrix
from sputils import getdtype, isshape, issequence, isscalarlike

class lil_matrix(spmatrix):
    """Row-based linked list sparse matrix

    This is an efficient structure for constructing sparse
    matrices incrementally.

    This can be instantiated in several ways:
        lil_matrix(D)
            with a dense matrix or rank-2 ndarray D

        lil_matrix(S)
            with another sparse matrix S (equivalent to S.tocsc())

        lil_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

    Notes
    -----

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

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        spmatrix.__init__(self)
        self.dtype = getdtype(dtype, arg1, default=float)

        # First get the shape
        if isspmatrix(arg1):
            if isspmatrix_lil(arg1) and copy:
                A = arg1.copy()
            else:
                A = arg1.tolil()

            if dtype is not None:
                A = A.astype(dtype)

            self.shape = A.shape
            self.dtype = A.dtype
            self.rows  = A.rows
            self.data  = A.data
        elif isinstance(arg1,tuple):
            if isshape(arg1):
                if shape is not None:
                    raise ValueError('invalid use of shape parameter')
                M, N = arg1
                self.shape = (M,N)
                self.rows = np.empty((M,), dtype=object)
                self.data = np.empty((M,), dtype=object)
                for i in range(M):
                    self.rows[i] = []
                    self.data[i] = []
            else:
                raise TypeError('unrecognized lil_matrix constructor usage')
        else:
            #assume A is dense
            try:
                A = np.asmatrix(arg1)
            except TypeError:
                raise TypeError('unsupported matrix type')
            else:
                from csr import csr_matrix
                A = csr_matrix(A, dtype=dtype).tolil()

                self.shape = A.shape
                self.dtype = A.dtype
                self.rows  = A.rows
                self.data  = A.data

    def __iadd__(self,other):
        self[:,:] = self + other
        return self

    def __isub__(self,other):
        self[:,:] = self - other
        return self

    def __imul__(self,other):
        if isscalarlike(other):
            self[:,:] = self * other
            return self
        else:
            raise NotImplementedError

    def __itruediv__(self,other):
        if isscalarlike(other):
            self[:,:] = self / other
            return self
        else:
            raise NotImplementedError

    # Whenever the dimensions change, empty lists should be created for each
    # row

    def getnnz(self):
        return sum([len(rowvals) for rowvals in self.data])
    nnz = property(fget=getnnz)

    def __str__(self):
        val = ''
        for i, row in enumerate(self.rows):
            for pos, j in enumerate(row):
                val += "  %s\t%s\n" % (str((i, j)), str(self.data[i][pos]))
        return val[:-1]

    def getrowview(self, i):
        """Returns a view of the 'i'th row (without copying).
        """
        new = lil_matrix((1, self.shape[1]), dtype=self.dtype)
        new.rows[0] = self.rows[i]
        new.data[0] = self.data[i]
        return new

    def getrow(self, i):
        """Returns a copy of the 'i'th row.
        """
        new = lil_matrix((1, self.shape[1]), dtype=self.dtype)
        new.rows[0] = self.rows[i][:]
        new.data[0] = self.data[i][:]
        return new

    def _get1(self, i, j):

        if i < 0:
            i += self.shape[0]
        if i < 0 or i >= self.shape[0]:
            raise IndexError('row index out of bounds')

        if j < 0:
            j += self.shape[1]
        if j < 0 or j >= self.shape[1]:
            raise IndexError('column index out of bounds')

        row  = self.rows[i]
        data = self.data[i]

        pos = bisect_left(row, j)
        if pos != len(data) and row[pos] == j:
            return data[pos]
        else:
            return 0

    def _slicetoseq(self, j, shape):
        if j.start is not None and j.start < 0:
            start =  shape + j.start
        elif j.start is None:
            start = 0
        else:
            start = j.start
        if j.stop is not None and j.stop < 0:
            stop = shape + j.stop
        elif j.stop is None:
            stop = shape
        else:
            stop = j.stop
        j = range(start, stop, j.step or 1)
        return j


    def __getitem__(self, index):
        """Return the element(s) index=(i, j), where j may be a slice.
        This always returns a copy for consistency, since slices into
        Python lists return copies.
        """
        try:
            i, j = index
        except (AssertionError, TypeError):
            raise IndexError('invalid index')

        if np.isscalar(i):
            if np.isscalar(j):
                return self._get1(i, j)
            if isinstance(j, slice):
                j = self._slicetoseq(j, self.shape[1])
            if issequence(j):
                return self.__class__([[self._get1(i, jj) for jj in j]])
        elif issequence(i) and issequence(j):
            return self.__class__([[self._get1(ii, jj) for (ii, jj) in zip(i, j)]])
        elif issequence(i) or isinstance(i, slice):
            if isinstance(i, slice):
                i = self._slicetoseq(i, self.shape[0])
            if np.isscalar(j):
                return self.__class__([[self._get1(ii, j)] for ii in i])
            if isinstance(j, slice):
                j = self._slicetoseq(j, self.shape[1])
            if issequence(j):
                return self.__class__([[self._get1(ii, jj) for jj in j] for ii in i])
        else:
            raise IndexError

    def _insertat2(self, row, data, j, x):
        """ helper for __setitem__: insert a value in the given row/data at
        column j. """

        if j < 0: #handle negative column indices
            j += self.shape[1]

        if j < 0 or j >= self.shape[1]:
            raise IndexError('column index out of bounds')
            
        if not np.isscalar(x):
            raise ValueError('setting an array element with a sequence')

        try:
            x = self.dtype.type(x)
        except:
            raise TypeError('Unable to convert value (%s) to dtype [%s]' % (x,self.dtype.name))

        pos = bisect_left(row, j)
        if x != 0:
            if pos == len(row):
                row.append(j)
                data.append(x)
            elif row[pos] != j:
                row.insert(pos, j)
                data.insert(pos, x)
            else:
                data[pos] = x
        else:
            if pos < len(row) and row[pos] == j:
                del row[pos]
                del data[pos]

    def _setitem_setrow(self, row, data, j, xrow, xdata, xcols):
        if isinstance(j, slice):
            j = self._slicetoseq(j, self.shape[1])
        if issequence(j):
            if xcols == len(j):
                for jj, xi in zip(j, xrange(xcols)):
                   pos = bisect_left(xrow, xi)
                   if pos != len(xdata) and xrow[pos] == xi:
                       self._insertat2(row, data, jj, xdata[pos])
                   else:
                       self._insertat2(row, data, jj, 0)
            elif xcols == 1:           # OK, broadcast across row
                if len(xdata) > 0 and xrow[0] == 0:
                    val = xdata[0]
                else:
                    val = 0
                for jj in j:
                    self._insertat2(row, data, jj,val)
            else:
                raise IndexError('invalid index')
        elif np.isscalar(j):
            if not xcols == 1:
                raise ValueError('array dimensions are not compatible for copy')
            if len(xdata) > 0 and xrow[0] == 0:
                self._insertat2(row, data, j, xdata[0])
            else:
                self._insertat2(row, data, j, 0)
        else:
            raise ValueError('invalid column value: %s' % str(j))

    def __setitem__(self, index, x):
        try:
            i, j = index
        except (ValueError, TypeError):
            raise IndexError('invalid index')

        # shortcut for common case of single entry assign:
        if np.isscalar(x) and np.isscalar(i) and np.isscalar(j):
            self._insertat2(self.rows[i], self.data[i], j, x)
            return

        # shortcut for common case of full matrix assign:
        if isspmatrix(x):
          if isinstance(i, slice) and i == slice(None) and \
             isinstance(j, slice) and j == slice(None):
               x = lil_matrix(x)
               self.rows = x.rows
               self.data = x.data
               return

        if isinstance(i, tuple):       # can't index lists with tuple
            i = list(i)

        if np.isscalar(i):
            rows = [self.rows[i]]
            datas = [self.data[i]]
        else:
            rows = self.rows[i]
            datas = self.data[i]

        x = lil_matrix(x, copy=False)
        xrows, xcols = x.shape
        if xrows == len(rows):    # normal rectangular copy
            for row, data, xrow, xdata in zip(rows, datas, x.rows, x.data):
                self._setitem_setrow(row, data, j, xrow, xdata, xcols)
        elif xrows == 1:          # OK, broadcast down column
            for row, data in zip(rows, datas):
                self._setitem_setrow(row, data, j, x.rows[0], x.data[0], xcols)

        # needed to pass 'test_lil_sequence_assignement' unit test:
        # -- set row from column of entries --
        elif xcols == len(rows):
            x = x.T
            for row, data, xrow, xdata in zip(rows, datas, x.rows, x.data):
                self._setitem_setrow(row, data, j, xrow, xdata, xrows)
        else:
            raise IndexError('invalid index')

    def _mul_scalar(self, other):
        if other == 0:
            # Multiply by zero: return the zero matrix
            new = lil_matrix(self.shape, dtype=self.dtype)
        else:
            new = self.copy()
            # Multiply this scalar by every element.
            new.data = np.array([[val*other for val in rowvals] for
                                  rowvals in new.data], dtype=object)
        return new

    def __truediv__(self, other):           # self / other
        if isscalarlike(other):
            new = self.copy()
            # Divide every element by this scalar
            new.data = np.array([[val/other for val in rowvals] for
                                  rowvals in new.data], dtype=object)
            return new
        else:
            return self.tocsr() / other

## This code doesn't work with complex matrices
#    def multiply(self, other):
#        """Point-wise multiplication by another lil_matrix.
#
#        """
#        if np.isscalar(other):
#            return self.__mul__(other)
#
#        if isspmatrix_lil(other):
#            reference,target = self,other
#
#            if reference.shape != target.shape:
#                raise ValueError("Dimensions do not match.")
#
#            if len(reference.data) > len(target.data):
#                reference,target = target,reference
#
#            new = lil_matrix(reference.shape)
#            for r,row in enumerate(reference.rows):
#                tr = target.rows[r]
#                td = target.data[r]
#                rd = reference.data[r]
#                L = len(tr)
#                for c,column in enumerate(row):
#                    ix = bisect_left(tr,column)
#                    if ix < L and tr[ix] == column:
#                        new.rows[r].append(column)
#                        new.data[r].append(rd[c] * td[ix])
#            return new
#        else:
#            raise ValueError("Point-wise multiplication only allowed "
#                             "with another lil_matrix.")

    def copy(self):
        from copy import deepcopy
        new = lil_matrix(self.shape, dtype=self.dtype)
        new.data = deepcopy(self.data)
        new.rows = deepcopy(self.rows)
        return new

    def reshape(self,shape):
        new = lil_matrix(shape, dtype=self.dtype)
        j_max = self.shape[1]
        for i,row in enumerate(self.rows):
            for col,j in enumerate(row):
                new_r,new_c = np.unravel_index(i*j_max + j,shape)
                new[new_r,new_c] = self[i,j]
        return new

    def toarray(self):
        d = np.zeros(self.shape, dtype=self.dtype)
        for i, row in enumerate(self.rows):
            for pos, j in enumerate(row):
                d[i, j] = self.data[i][pos]
        return d

    def transpose(self):
        return self.tocsr().transpose().tolil()

    def tolil(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self):
        """ Return Compressed Sparse Row format arrays for this matrix.
        """

        indptr = np.asarray([len(x) for x in self.rows], dtype=np.intc)
        indptr = np.concatenate( (np.array([0], dtype=np.intc), np.cumsum(indptr)) )

        nnz = indptr[-1]

        indices = []
        for x in self.rows:
            indices.extend(x)
        indices = np.asarray(indices, dtype=np.intc)

        data = []
        for x in self.data:
            data.extend(x)
        data = np.asarray(data, dtype=self.dtype)

        from csr import csr_matrix
        return csr_matrix((data, indices, indptr), shape=self.shape)

    def tocsc(self):
        """ Return Compressed Sparse Column format arrays for this matrix.
        """
        return self.tocsr().tocsc()


from sputils import _isinstance

def isspmatrix_lil( x ):
    return _isinstance(x, lil_matrix)
