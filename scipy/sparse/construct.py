""" Functions to construct sparse matrices
"""


__all__ = [ 'spdiags','speye','spidentity','spkron', 'lil_eye', 'lil_diags' ]

import itertools

import numpy
from numpy import ones, clip, array, arange, intc

from sparse import csr_matrix, csc_matrix, coo_matrix, \
        dok_matrix, lil_matrix
from sparse import isspmatrix, isspmatrix_csr, isspmatrix_csc
import sparsetools

def _spdiags_tosub(diag_num, a, b):
    part1 = where(less(diag_num, a), abs(diag_num-a), 0)
    part2 = where(greater(diag_num, b), abs(diag_num-b), 0)
    return part1+part2

# Note: sparsetools only offers diagonal -> CSC matrix conversion functions,
# not to CSR
def spdiags(diags, offsets, M, N):
    """Return a sparse matrix in CSC format given its diagonals.

    B = spdiags(diags, offsets, M, N)

    Inputs:
        diags  --  rows contain diagonal values
        offsets -- diagonals to set (0 is main)
        M, N    -- sparse matrix returned is M X N
    """
    #    diags = array(transpose(diags), copy=True)
    diags = array(diags, copy = True)
    if diags.dtype.char not in 'fdFD':
        diags = diags.astype('d')
    if not hasattr(offsets, '__len__' ):
        offsets = (offsets,)
    offsets = array(offsets, copy=False, dtype=intc)
    assert(len(offsets) == diags.shape[0])
    indptr, rowind, data = sparsetools.spdiags(M, N, len(offsets), offsets, diags)
    return csc_matrix((data, rowind, indptr), (M, N))



def spidentity(n, dtype='d'):
    """
    spidentity( n ) returns the identity matrix of shape (n, n) stored
    in CSC sparse matrix format.
    """
    return csc_matrix((ones(n,dtype=dtype),arange(n),arange(n+1)),(n,n))


def speye(n, m, k = 0, dtype = 'd'):
    """
    speye(n, m) returns a (n x m) matrix stored
    in CSC sparse matrix format, where the  k-th diagonal is all ones,
    and everything else is zeros.
    """
    diags = ones((1, n), dtype = dtype)
    return spdiags(diags, k, n, m)

def spkron(a,b):
    """kronecker product of sparse matrices a and b

    *Parameters*:
        a,b : sparse matrices
            E.g. csr_matrix, csc_matrix, coo_matrix, etc.

    *Returns*:
        coo_matrix
            kronecker product in COOrdinate format

    *Example*:
    -------

    >>> a = csr_matrix(array([[0,2],[5,0]]))
    >>> b = csr_matrix(array([[1,2],[3,4]]))
    >>> spkron(a,b).todense()
    matrix([[  0.,   0.,   2.,   4.],
            [  0.,   0.,   6.,   8.],
            [  5.,  10.,   0.,   0.],
            [ 15.,  20.,   0.,   0.]])

    """
    if not isspmatrix(a) and isspmatrix(b):
        raise ValueError,'expected sparse matrix'

    a,b = a.tocoo(),b.tocoo()
    output_shape = (a.shape[0]*b.shape[0],a.shape[1]*b.shape[1])

    if a.nnz == 0 or b.nnz == 0:
        # kronecker product is the zero matrix
        return coo_matrix( output_shape )


    # expand entries of a into blocks
    row  = a.row.repeat(b.nnz)
    col  = a.col.repeat(b.nnz)
    data = a.data.repeat(b.nnz)

    row *= b.shape[0]
    col *= b.shape[1]

    # increment block indices
    row,col = row.reshape(-1,b.nnz),col.reshape(-1,b.nnz)
    row += b.row
    col += b.col
    row,col = row.reshape(-1),col.reshape(-1)

    # compute block entries
    data = data.reshape(-1,b.nnz) * b.data
    data = data.reshape(-1)

    return coo_matrix((data,(row,col)), dims=output_shape)




def lil_eye((r,c), k=0, dtype=float):
    """Generate a lil_matrix of dimensions (r,c) with the k-th
    diagonal set to 1.

    :Parameters:
        r,c : int
            Row and column-dimensions of the output.
        k : int
            Diagonal offset.  In the output matrix,
            out[m,m+k] == 1 for all m.
        dtype : dtype
            Data-type of the output array.

    """
    out = lil_matrix((r,c),dtype=dtype)
    for c in xrange(clip(k,0,c),clip(r+k,0,c)):
        out.rows[c-k].append(c)
        out.data[c-k].append(1)
    return out

def lil_diags(diags,offsets,(m,n),dtype=float):
    """Generate a lil_matrix with the given diagonals.

    :Parameters:
        diags : list of list of values e.g. [[1,2,3],[4,5]]
            Values to be placed on each indicated diagonal.
        offsets : list of ints
            Diagonal offsets.  This indicates the diagonal on which
            the given values should be placed.
        (r,c) : tuple of ints
            Row and column dimensions of the output.
        dtype : dtype
           Output data-type.

    Example:
    -------

    >>> lil_diags([[1,2,3],[4,5],[6]],[0,1,2],(3,3)).todense()
    matrix([[ 1.,  4.,  6.],
            [ 0.,  2.,  5.],
            [ 0.,  0.,  3.]])

    """
    offsets_unsorted = list(offsets)
    diags_unsorted = list(diags)
    if len(diags) != len(offsets):
        raise ValueError("Number of diagonals provided should "
                         "agree with offsets.")

    sort_indices = numpy.argsort(offsets_unsorted)
    diags = [diags_unsorted[k] for k in sort_indices]
    offsets = [offsets_unsorted[k] for k in sort_indices]

    for i,k in enumerate(offsets):
        if len(diags[i]) < m-abs(k):
            raise ValueError("Not enough values specified to fill "
                             "diagonal %s." % k)

    out = lil_matrix((m,n),dtype=dtype)
    for k,diag in itertools.izip(offsets,diags):
        for ix,c in enumerate(xrange(clip(k,0,n),clip(m+k,0,n))):
            out.rows[c-k].append(c)
            out.data[c-k].append(diag[ix])
    return out


