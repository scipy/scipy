"""LInked List sparse matrix class
"""

__docformat__ = "restructuredtext en"

__all__ = ['lil_matrix','isspmatrix_lil']

import copy
from bisect import bisect_left

import numpy
from numpy import isscalar, asmatrix, asarray, intc, concatenate, array, \
        cumsum, zeros, unravel_index

from base import spmatrix, isspmatrix
from sputils import getdtype,isshape,issequence,isscalarlike

class lil_matrix(spmatrix):
    """Row-based linked list matrix


    This can be instantiated in several ways:
        csc_matrix(D)
            with a dense matrix or rank-2 ndarray D

        csc_matrix(S)
            with another sparse matrix S (equivalent to S.tocsc())

        csc_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        csc_matrix((data, ij), [shape=(M, N)])
            where ``data`` and ``ij`` satisfy ``a[ij[0, k], ij[1, k]] = data[k]``

    Notes
    -----

    Advantages of the LIL format
        - supports flexible slicing
        - changes to the matrix sparsity structure are efficient

    Disadvantages of the LIL format
        - arithmetic operations LIL + LIL are slow (consider CSR or CSC)
        - slow column slicing (consider CSC)
        - matrix vector products are slower than CSR/CSC

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

    def __init__(self, A=None, shape=None, dtype=None, copy=False):
        """ Create a new list-of-lists sparse matrix.  An optional
        argument A is accepted, which initializes the lil_matrix with it.
        This can be a tuple of dimensions (M, N) or a dense array /
        matrix to copy, or a sparse matrix.
        """
        spmatrix.__init__(self)
        self.dtype = getdtype(dtype, A, default=float)

        # First get the shape
        if A is None:
            if not isshape(shape):
                raise TypeError("need a valid shape")
            M, N = shape
            self.shape = (M,N)
            self.rows = numpy.empty((M,), dtype=object)
            self.data = numpy.empty((M,), dtype=object)
            for i in range(M):
                self.rows[i] = []
                self.data[i] = []
        elif isspmatrix(A):
            if isspmatrix_lil(A) and copy:
                A = A.copy()
            else:
                A = A.tolil()
            self.shape = A.shape
            self.dtype = A.dtype
            self.rows  = A.rows
            self.data  = A.data
        elif isinstance(A,tuple):
            if isshape(A):
                if shape is not None:
                    raise ValueError('invalid use of shape parameter')
                M, N = A
                self.shape = (M,N)
                self.rows = numpy.empty((M,), dtype=object)
                self.data = numpy.empty((M,), dtype=object)
                for i in range(M):
                    self.rows[i] = []
                    self.data[i] = []
            else:
                raise TypeError,'unrecognized lil_matrix constructor usage'
        else:
            #assume A is dense
            try:
                A = asmatrix(A)
            except TypeError:
                raise TypeError, "unsupported matrix type"
            else:
                from csr import csr_matrix
                A = csr_matrix(A).tolil()

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
        row = self.rows[i]
        data = self.data[i]

        if j < 0:
            j += self.shape[1]

        if j < 0 or j > self.shape[1]:
            raise IndexError,'column index out of bounds'

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
            raise IndexError, "invalid index"

        if isscalar(i):
            if isscalar(j):
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
            if isscalar(j):
                return self.__class__([[self._get1(ii, j)] for ii in i])
            if isinstance(j, slice):
                j = self._slicetoseq(j, self.shape[1])
            if issequence(j):
                return self.__class__([[self._get1(ii, jj) for jj in j] for ii in i])
        else:
            raise IndexError


    def _insertat(self, i, j, x):
        """ helper for __setitem__: insert a value at (i,j) where i, j and x
        are all scalars """
        row = self.rows[i]
        data = self.data[i]
        self._insertat2(row, data, j, x)

    def _insertat2(self, row, data, j, x):
        """ helper for __setitem__: insert a value in the given row/data at
        column j. """

        if j < 0: #handle negative column indices
            j += self.shape[1]

        if j < 0 or j >= self.shape[1]:
            raise IndexError,'column index out of bounds'

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


    def _insertat3(self, row, data, j, x):
        """ helper for __setitem__ """
        if isinstance(j, slice):
            j = self._slicetoseq(j, self.shape[1])
        if issequence(j):
            if isinstance(x, spmatrix):
                x = x.todense()
            x = numpy.asarray(x).squeeze()
            if isscalar(x) or x.size == 1:
                for jj in j:
                    self._insertat2(row, data, jj, x)
            else:
                # x must be one D. maybe check these things out
                for jj, xx in zip(j, x):
                    self._insertat2(row, data, jj, xx)
        elif isscalar(j):
            self._insertat2(row, data, j, x)
        else:
            raise ValueError, "invalid column value: %s" % str(j)


    def __setitem__(self, index, x):
        if isscalar(x):
            x = self.dtype.type(x)
        elif not isinstance(x, spmatrix):
            x = numpy.asarray(x, dtype=self.dtype)

        try:
            i, j = index
        except (ValueError, TypeError):
            raise IndexError, "invalid index"

        if isscalar(i):
            row = self.rows[i]
            data = self.data[i]
            self._insertat3(row, data, j, x)
        elif issequence(i) and issequence(j):
            if isscalar(x):
                for ii, jj in zip(i, j):
                    self._insertat(ii, jj, x)
            else:
                for ii, jj, xx in zip(i, j, x):
                    self._insertat(ii, jj, xx)
        elif isinstance(i, slice) or issequence(i):
            rows = self.rows[i]
            datas = self.data[i]
            if isscalar(x):
                for row, data in zip(rows, datas):
                    self._insertat3(row, data, j, x)
            else:
                for row, data, xx in zip(rows, datas, x):
                    self._insertat3(row, data, j, xx)
        else:
            raise ValueError, "invalid index value: %s" % str((i, j))

    def __mul__(self, other):           # self * other
        if isscalarlike(other):
            if other == 0:
                # Multiply by zero: return the zero matrix
                new = lil_matrix(shape=self.shape, dtype=self.dtype)
            else:
                new = self.copy()
                # Multiply this scalar by every element.
                new.data = numpy.array([[val*other for val in rowvals] for
                                        rowvals in new.data], dtype=object)
            return new
        else:
            return self.dot(other)

    def __truediv__(self, other):           # self / other
        if isscalarlike(other):
            new = self.copy()
            # Divide every element by this scalar
            new.data = numpy.array([[val/other for val in rowvals] for
                                    rowvals in new.data], dtype=object)
            return new
        else:
            return self.tocsr() / other

    def multiply(self, other):
        """Point-wise multiplication by another lil_matrix.

        """
        if isscalar(other):
            return self.__mul__(other)

        if isspmatrix_lil(other):
            reference,target = self,other

            if reference.shape != target.shape:
                raise ValueError("Dimensions do not match.")

            if len(reference.data) > len(target.data):
                reference,target = target,reference

            new = lil_matrix(reference.shape)
            for r,row in enumerate(reference.rows):
                tr = target.rows[r]
                td = target.data[r]
                rd = reference.data[r]
                L = len(tr)
                for c,column in enumerate(row):
                    ix = bisect_left(tr,column)
                    if ix < L and tr[ix] == column:
                        new.rows[r].append(column)
                        new.data[r].append(rd[c] * td[ix])
            return new
        else:
            raise ValueError("Point-wise multiplication only allowed "
                             "with another lil_matrix.")

    def copy(self):
        new = lil_matrix(self.shape, dtype=self.dtype)
        new.data = copy.deepcopy(self.data)
        new.rows = copy.deepcopy(self.rows)
        return new

    def reshape(self,shape):
        new = lil_matrix(shape,dtype=self.dtype)
        j_max = self.shape[1]
        for i,row in enumerate(self.rows):
            for col,j in enumerate(row):
                new_r,new_c = unravel_index(i*j_max + j,shape)
                new[new_r,new_c] = self[i,j]
        return new

    def __add__(self, other):
        if isscalar(other) and other != 0:
            raise ValueError("Refusing to destroy sparsity. "
                             "Use x.todense() + c instead.")
        else:
            return spmatrix.__add__(self, other)

    def __rmul__(self, other):          # other * self
        if isscalarlike(other):
            # Multiplication by a scalar is symmetric
            return self.__mul__(other)
        else:
            return spmatrix.__rmul__(self, other)


    def toarray(self):
        d = zeros(self.shape, dtype=self.dtype)
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

        indptr = asarray([len(x) for x in self.rows], dtype=intc)
        indptr = concatenate( ( array([0],dtype=intc), cumsum(indptr) ) )

        nnz = indptr[-1]

        indices = []
        for x in self.rows:
            indices.extend(x)
        indices = asarray(indices,dtype=intc)

        data = []
        for x in self.data:
            data.extend(x)
        data = asarray(data,dtype=self.dtype)

        from csr import csr_matrix
        return csr_matrix((data, indices, indptr), shape=self.shape)

    def tocsc(self):
        """ Return Compressed Sparse Column format arrays for this matrix.
        """
        return self.tocsr().tocsc()


from sputils import _isinstance

def isspmatrix_lil( x ):
    return _isinstance(x, lil_matrix)
