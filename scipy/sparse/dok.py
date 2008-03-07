"""Dictionary Of Keys based matrix"""

__docformat__ = "restructuredtext en"

__all__ = ['dok_matrix', 'isspmatrix_dok']

import operator
from itertools import izip

from numpy import asarray, asmatrix, intc, isscalar, array, matrix

from base import spmatrix,isspmatrix
from sputils import isdense, getdtype, isshape, isintlike, isscalarlike

class dok_matrix(spmatrix, dict):
    """Dictionary Of Keys based matrix.  This is an efficient
    structure for constructing sparse matrices 
    """
    def __init__(self, A=None, shape=None, dtype=None, copy=False):
        """ Create a new dictionary-of-keys sparse matrix.  An optional
        argument A is accepted, which initializes the dok_matrix with it.
        This can be a tuple of dimensions (M, N) or a (dense) array
        to copy.
        """
        #TODO deprecate argument A in favor of arg1 style

        dict.__init__(self)
        spmatrix.__init__(self)
        self.dtype = getdtype(dtype, A, default=float)
        if A is not None:
            if isinstance(A, tuple):
                # Interpret as dimensions
                if not isshape(A):
                    raise TypeError, "dimensions must be a 2-tuple of positive"\
                            " integers"
                self.shape = A
            elif isspmatrix(A):
                if isspmatrix_dok(A) and copy:
                    A = A.copy()
                else:
                    A = A.todok()
                self.update( A )
                self.shape = A.shape
                self.dtype = A.dtype
            else:
                #must be dense, convert to COO first, then to DOK
                try:
                    A = asarray(A)
                except:
                    raise ValueError, "unrecognized form for" \
                            " %s_matrix constructor" % self.format
                from coo import coo_matrix
                self.update( coo_matrix(A).todok() )
                self.shape = A.shape
                self.dtype = A.dtype

    def getnnz(self):
        return dict.__len__(self)
    nnz = property(fget=getnnz)

    def __len__(self):
        return dict.__len__(self)

    def get(self, key, default=0.):
        """This overrides the dict.get method, providing type checking
        but otherwise equivalent functionality.
        """
        try:
            i, j = key
            assert isintlike(i) and isintlike(j)
        except (AssertionError, TypeError, ValueError):
            raise IndexError, "index must be a pair of integers"
        try:
            assert not (i < 0 or i >= self.shape[0] or j < 0 or
                     j >= self.shape[1])
        except AssertionError:
            raise IndexError, "index out of bounds"
        return dict.get(self, key, default)

    def  __getitem__(self, key):
        """If key=(i,j) is a pair of integers, return the corresponding
        element.  If either i or j is a slice or sequence, return a new sparse
        matrix with just these elements.
        """
        try:
            i, j = key
        except (ValueError, TypeError):
            raise TypeError, "index must be a pair of integers or slices"


        # Bounds checking
        if isintlike(i):
            if i < 0:
                i += self.shape[0]
            if i < 0 or i >= self.shape[0]:
                raise IndexError, "index out of bounds"
        if isintlike(j):
            if j < 0:
                j += self.shape[1]
            if j < 0 or j >= self.shape[1]:
                raise IndexError, "index out of bounds"

        # First deal with the case where both i and j are integers
        if isintlike(i) and isintlike(j):
            return dict.get(self, key, 0.)
        else:
            # Either i or j is a slice, sequence, or invalid.  If i is a slice
            # or sequence, unfold it first and call __getitem__ recursively.

            if isinstance(i, slice):
                # Is there an easier way to do this?
                seq = xrange(i.start or 0, i.stop or self.shape[0], i.step or 1)
            elif operator.isSequenceType(i):
                seq = i
            else:
                # Make sure i is an integer. (But allow it to be a subclass of int).
                if not isintlike(i):
                    raise TypeError, "index must be a pair of integers or slices"
                seq = None
            if seq is not None:
                # i is a seq
                if isintlike(j):
                    # Create a new matrix of the correct dimensions
                    first = seq[0]
                    last = seq[-1]
                    if first < 0 or first >= self.shape[0] or last < 0 \
                                 or last >= self.shape[0]:
                        raise IndexError, "index out of bounds"
                    newshape = (last-first+1, 1)
                    new = dok_matrix(newshape)
                    # ** This uses linear time in the size m of dimension 0:
                    # new[0:seq[-1]-seq[0]+1, 0] = \
                    #         [self.get((element, j), 0) for element in seq]
                    # ** Instead just add the non-zero elements.  This uses
                    # ** linear time in the number of non-zeros:
                    for (ii, jj) in self.keys():
                        if jj == j and ii >= first and ii <= last:
                            dict.__setitem__(new, (ii-first, 0), \
                                             dict.__getitem__(self, (ii,jj)))
                else:
                    ###################################
                    # We should reshape the new matrix here!
                    ###################################
                    raise NotImplementedError, "fancy indexing supported over" \
                            " one axis only"
                return new

            # Below here, j is a sequence, but i is an integer
            if isinstance(j, slice):
                # Is there an easier way to do this?
                seq = xrange(j.start or 0, j.stop or self.shape[1], j.step or 1)
            elif operator.isSequenceType(j):
                seq = j
            else:
                # j is not an integer
                raise TypeError, "index must be a pair of integers or slices"

            # Create a new matrix of the correct dimensions
            first = seq[0]
            last = seq[-1]
            if first < 0 or first >= self.shape[1] or last < 0 \
                         or last >= self.shape[1]:
                raise IndexError, "index out of bounds"
            newshape = (1, last-first+1)
            new = dok_matrix(newshape)
            # ** This uses linear time in the size n of dimension 1:
            # new[0, 0:seq[-1]-seq[0]+1] = \
            #         [self.get((i, element), 0) for element in seq]
            # ** Instead loop over the non-zero elements.  This is slower
            # ** if there are many non-zeros
            for (ii, jj) in self.keys():
                if ii == i and jj >= first and jj <= last:
                    dict.__setitem__(new, (0, jj-first), \
                                     dict.__getitem__(self, (ii,jj)))
            return new


    def __setitem__(self, key, value):
        try:
            assert len(key) == 2
        except (AssertionError, TypeError):
            raise TypeError, "index must be a pair of integers, slices, or" \
                    " sequences"
        i, j = key


        # First deal with the case where both i and j are integers
        if isintlike(i) and isintlike(j):
            if i < 0:
                i += self.shape[0]
            if j < 0:
                j += self.shape[1]

            if i < 0 or i >= self.shape[0] or j < 0 or j >= self.shape[1]:
                raise IndexError, "index out of bounds"
            if isintlike(value) and value == 0:
                if key in self.keys():  # get rid of it something already there
                    del self[key]
            else:
                # Ensure value is a single element, not a sequence
                if isinstance(value, float) or isintlike(value) or \
                        isinstance(value, complex):
                    dict.__setitem__(self, (i,j), self.dtype.type(value))
                    newrows = max(self.shape[0], int(key[0])+1)
                    newcols = max(self.shape[1], int(key[1])+1)
                    self.shape = (newrows, newcols)
                else:
                    raise TypeError, "cannot set matrix element to non-scalar"
            return                 # done
        else:
            # Either i or j is a slice, sequence, or invalid.  If i is a slice
            # or sequence, unfold it first and call __setitem__ recursively.
            if isinstance(i, slice):
                # Is there an easier way to do this?
                seq = xrange(i.start or 0, i.stop or self.shape[0], i.step or 1)
            elif operator.isSequenceType(i):
                seq = i
            else:
                # Make sure i is an integer. (But allow it to be a subclass of int).
                if not isintlike(i):
                    raise TypeError, "index must be a pair of integers or slices"
                seq = None
            if seq is not None:
                # First see if 'value' is another dok_matrix of the appropriate
                # dimensions
                if isinstance(value, dok_matrix):
                    if value.shape[1] == 1:
                        for element in seq:
                            self[element, j] = value[element, 0]
                    else:
                        raise NotImplementedError, "setting a 2-d slice of" \
                                " a dok_matrix is not yet supported"
                elif isscalar(value):
                    for element in seq:
                        self[element, j] = value
                else:
                    # See if value is a sequence
                    try:
                        if len(seq) != len(value):
                            raise ValueError, "index and value ranges must" \
                                              " have the same length"
                    except TypeError:
                        # Not a sequence
                        raise TypeError, "unsupported type for" \
                                         " dok_matrix.__setitem__"

                    # Value is a sequence
                    for element, val in izip(seq, value):
                        self[element, j] = val   # don't use dict.__setitem__
                            # here, since we still want to be able to delete
                            # 0-valued keys, do type checking on 'val' (e.g. if
                            # it's a rank-1 dense array), etc.
            else:
                # Process j
                if isinstance(j, slice):
                    seq = xrange(j.start or 0, j.stop or self.shape[1], j.step or 1)
                elif operator.isSequenceType(j):
                    seq = j
                else:
                    # j is not an integer
                    raise TypeError, "index must be a pair of integers or slices"

                # First see if 'value' is another dok_matrix of the appropriate
                # dimensions
                if isinstance(value, dok_matrix):
                    if value.shape[0] == 1:
                        for element in seq:
                            self[i, element] = value[0, element]
                    else:
                        raise NotImplementedError, "setting a 2-d slice of" \
                                " a dok_matrix is not yet supported"
                elif isscalar(value):
                    for element in seq:
                        self[i, element] = value
                else:
                    # See if value is a sequence
                    try:
                        if len(seq) != len(value):
                            raise ValueError, "index and value ranges must have" \
                                              " the same length"
                    except TypeError:
                        # Not a sequence
                        raise TypeError, "unsupported type for dok_matrix.__setitem__"
                    else:
                        for element, val in izip(seq, value):
                            self[i, element] = val


    def __add__(self, other):
        # First check if argument is a scalar
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Add this scalar to every element.
            M, N = self.shape
            for i in xrange(M):
                for j in xrange(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
            #new.dtype.char = self.dtype.char
        elif isinstance(other, dok_matrix):
            if other.shape != self.shape:
                raise ValueError, "matrix dimensions are not equal"
            # We could alternatively set the dimensions to the the largest of
            # the two matrices to be summed.  Would this be a good idea?
            new = dok_matrix(self.shape, dtype=self.dtype)
            new.update(self)
            for key in other.keys():
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = self.todense() + other
        else:
            raise TypeError, "data type not understood"
        return new

    def __radd__(self, other):
        # First check if argument is a scalar
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Add this scalar to every element.
            M, N = self.shape
            for i in xrange(M):
                for j in xrange(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
        elif isinstance(other, dok_matrix):
            if other.shape != self.shape:
                raise ValueError, "matrix dimensions are not equal"
            new = dok_matrix(self.shape, dtype=self.dtype)
            new.update(self)
            for key in other:
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = other + self.todense()
        else:
            raise TypeError, "data type not understood"
        return new

    def __neg__(self):
        new = dok_matrix(self.shape, dtype=self.dtype)
        for key in self.keys():
            new[key] = -self[key]
        return new

    def __mul__(self, other):           # self * other
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                new[key] = val * other
            #new.dtype.char = self.dtype.char
            return new
        else:
            return self.dot(other)

    def __imul__(self, other):           # self * other
        if isscalarlike(other):
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                self[key] = val * other
            #new.dtype.char = self.dtype.char
            return self
        else:
            return NotImplementedError

    def __rmul__(self, other):          # other * self
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                new[key] = other * val
            #new.dtype.char = self.dtype.char
            return new
        else:
            # Don't use asarray unless we have to
            try:
                tr = other.transpose()
            except AttributeError:
                tr = asarray(other).transpose()
            return self.transpose().dot(tr).transpose()

    def __truediv__(self, other):           # self * other
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                new[key] = val / other
            #new.dtype.char = self.dtype.char
            return new
        else:
            return self.tocsr() / other


    def __itruediv__(self, other):           # self * other
        if isscalarlike(other):
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                self[key] = val / other
            return self
        else:
            return NotImplementedError

    # What should len(sparse) return? For consistency with dense matrices,
    # perhaps it should be the number of rows?  For now it returns the number
    # of non-zeros.

    def transpose(self):
        """ Return the transpose
        """
        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype)
        for key, value in self.iteritems():
            new[key[1], key[0]] = value
        return new

    def conjtransp(self):
        """ Return the conjugate transpose
        """
        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype)
        for key, value in self.iteritems():
            new[key[1], key[0]] = conj(value)
        return new

    def copy(self):
        new = dok_matrix(self.shape, dtype=self.dtype)
        new.update(self)
        new.shape = self.shape
        return new

    def take(self, cols_or_rows, columns=1):
        # Extract columns or rows as indictated from matrix
        # assume cols_or_rows is sorted
        new = dok_matrix(dtype=self.dtype)    # what should the dimensions be ?!
        indx = int((columns == 1))
        N = len(cols_or_rows)
        if indx: # columns
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[1])
                if num < N:
                    newkey = (key[0], num)
                    new[newkey] = self[key]
        else:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[0])
                if num < N:
                    newkey = (num, key[1])
                    new[newkey] = self[key]
        return new

    def split(self, cols_or_rows, columns=1):
        # Similar to take but returns two arrays, the extracted columns plus
        # the resulting array.  Assumes cols_or_rows is sorted
        base = dok_matrix()
        ext = dok_matrix()
        indx = int((columns == 1))
        if indx:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[1])
                if cols_or_rows[num] == key[1]:
                    newkey = (key[0], num)
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0], key[1]-num)
                    base[newkey] = self[key]
        else:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[0])
                if cols_or_rows[num] == key[0]:
                    newkey = (num, key[1])
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0]-num, key[1])
                    base[newkey] = self[key]
        return base, ext


    def matvec(self, other):
        if isdense(other):
            if other.shape[0] != self.shape[1]:
                raise ValueError, "dimensions do not match"
            new = [0] * self.shape[0]
            for key in self.keys():
                new[int(key[0])] += self[key] * other[int(key[1]), ...]
            new = array(new)
            if isinstance(other, matrix):
                new = asmatrix(new)
                # Do we need to return the transpose?
                if other.shape[1] == 1:
                    new = new.T
            return new
        else:
            raise TypeError, "need a dense vector"

    def rmatvec(self, other, conjugate=True):
        if isdense(other):
            if other.shape[-1] != self.shape[0]:
                raise ValueError, "dimensions do not match"
            new = [0] * self.shape[1]
            for key in self.keys():
                new[int(key[1])] += other[..., int(key[0])] * conj(self[key])
            new = array(new)
            if isinstance(other, matrix):
                new = asmatrix(new)
                # Do we need to return the transpose?
                if other.shape[1] == 1:
                    new = new.T
            return new
        else:
            raise TypeError, "need a dense vector"

    def tocoo(self):
        """ Return a copy of this matrix in COOrdinate format"""
        from coo import coo_matrix
        if self.nnz == 0:
            return coo_matrix(self.shape, dtype=self.dtype)
        else:
            data    = asarray(self.values(), dtype=self.dtype)
            indices = asarray(self.keys(), dtype=intc).T
            return coo_matrix((data,indices),shape=self.shape,dtype=self.dtype)

    def todok(self,copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self):
        """ Return a copy of this matrix in Compressed Sparse Row format"""
        return self.tocoo().tocsr()

    def tocsc(self):
        """ Return a copy of this matrix in Compressed Sparse Column format"""
        return self.tocoo().tocsc()

    def toarray(self):
        return self.tocoo().toarray()

    def resize(self, shape):
        """ Resize the matrix to dimensions given by 'shape', removing any
        non-zero elements that lie outside.
        """
        if not isshape(shape):
            raise TypeError, "dimensions must be a 2-tuple of positive"\
                             " integers"
        newM, newN = shape
        M, N = self.shape
        if newM < M or newN < N:
            # Remove all elements outside new dimensions
            for (i, j) in self.keys():
                if i >= newM or j >= newN:
                    del self[i, j]
        self.shape = shape



from sputils import _isinstance

def isspmatrix_dok(x):
    return _isinstance(x, dok_matrix)

