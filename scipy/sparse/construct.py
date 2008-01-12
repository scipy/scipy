""" Functions to construct sparse matrices
"""


__all__ = [ 'spdiags','speye','spidentity','spkron', 'lil_eye', 'lil_diags' ]

from itertools import izip
from warnings import warn

import numpy
from numpy import ones, clip, array, arange, intc

from csr import csr_matrix, isspmatrix_csr
from csc import csc_matrix, isspmatrix_csc
from bsr import bsr_matrix
from coo import coo_matrix
from dok import dok_matrix
from lil import lil_matrix
from dia import dia_matrix
from base import isspmatrix


def spdiags(data, diags, m, n, format=None):
    """Return a sparse matrix given its diagonals.

    B = spdiags(diags, offsets, m, n)

    Parameters
    ==========
        - data   : matrix whose rows contain the diagonal values
        - diags  : diagonals to set 
            - k = 0 - the main diagonal
            - k > 0 - the k-th upper diagonal
            - k < 0 - the k-th lower diagonal
        - m, n   : dimensions of the result
        - format : format of the result (e.g. "csr")
            -  By default (format=None) an appropriate sparse matrix 
               format is returned.  This choice is subject to change.

    See Also
    ========
        The dia_matrix class which implements the DIAgonal format.

    Example
    =======
    >>> data = array([[1,2,3,4]]).repeat(3,axis=0)
    >>> diags = array([0,-1,2])
    >>> spdiags(data,diags,4,4).todense()
    matrix([[1, 0, 3, 0],
            [1, 2, 0, 4],
            [0, 2, 3, 0],
            [0, 0, 3, 4]])

    """
    return dia_matrix((data, diags), shape=(m,n)).asformat(format)

def spidentity(n, dtype='d', format=None):
    """spidentity(n) returns an (n x n) identity matrix"""
    if format in ['csr','csc']:
        indptr  = arange(n+1, dtype=intc)
        indices = arange(n, dtype=intc)
        data    = ones(n, dtype=dtype)
        cls = eval('%s_matrix' % format)
        return cls((data,indices,indptr),(n,n))
    elif format == 'coo':
        row  = arange(n, dtype=intc)
        col  = arange(n, dtype=intc)
        data = ones(n, dtype=dtype)
        cls = eval('%s_matrix' % format)
        return coo_matrix((data,(row,col)),(n,n))
    else:
        return spidentity( n, dtype=dtype, format='csr').asformat(format)

def speye(m, n, k=0, dtype='d', format=None):
    """speye(m, n) returns an (m x n) matrix where the k-th diagonal 
    is all ones and everything else is zeros.
    """
    diags = ones((1, m), dtype = dtype)
    return spdiags(diags, k, m, n).asformat(format)

def spkron(A, B, format=None):
    """kronecker product of sparse matrices A and B

    Parameters
    ==========
        A,B    : dense or sparse matrices
        format : format of the result (e.g. "csr")
            -  By default (format=None) an appropriate sparse matrix 
               format is returned.  This choice is subject to change.

    Returns
    =======
        kronecker product in a sparse matrix format

    Examples
    ========

    >>> A = csr_matrix(array([[0,2],[5,0]]))
    >>> B = csr_matrix(array([[1,2],[3,4]]))
    >>> spkron(A,B).todense()
    matrix([[  0.,   0.,   2.,   4.],
            [  0.,   0.,   6.,   8.],
            [  5.,  10.,   0.,   0.],
            [ 15.,  20.,   0.,   0.]])

    >>> spkron(A,[[1,2],[3,4]]).todense()
    matrix([[  0.,   0.,   2.,   4.],
            [  0.,   0.,   6.,   8.],
            [  5.,  10.,   0.,   0.],
            [ 15.,  20.,   0.,   0.]])

    """
    #TODO optimize for small dense B and CSR A -> BSR
    B = coo_matrix(B)

    
    if (format is None or format == "bsr") and 2*B.nnz >= B.shape[0] * B.shape[1]:
        #B is fairly dense, use BSR
        A = csr_matrix(A,copy=True)
        
        output_shape = (A.shape[0]*B.shape[0],A.shape[1]*B.shape[1])

        if A.nnz == 0 or B.nnz == 0:
            # kronecker product is the zero matrix
            return coo_matrix( output_shape )
        
        B = B.toarray()
        data = A.data.repeat(B.size).reshape(-1,B.shape[0],B.shape[1])
        data = data * B
        
        return bsr_matrix((data,A.indices,A.indptr),shape=output_shape)
    else:
        #use COO
        A = coo_matrix(A)
        output_shape = (A.shape[0]*B.shape[0],A.shape[1]*B.shape[1])

        if A.nnz == 0 or B.nnz == 0:
            # kronecker product is the zero matrix
            return coo_matrix( output_shape )

        # expand entries of a into blocks
        row  = A.row.repeat(B.nnz)
        col  = A.col.repeat(B.nnz)
        data = A.data.repeat(B.nnz)

        row *= B.shape[0]
        col *= B.shape[1]

        # increment block indices
        row,col = row.reshape(-1,B.nnz),col.reshape(-1,B.nnz)
        row += B.row
        col += B.col
        row,col = row.reshape(-1),col.reshape(-1)

        # compute block entries
        data = data.reshape(-1,B.nnz) * B.data
        data = data.reshape(-1)

        return coo_matrix((data,(row,col)), shape=output_shape).asformat(format)




def lil_eye((r,c), k=0, dtype='d'):
    """Generate a lil_matrix of dimensions (r,c) with the k-th
    diagonal set to 1.

    Parameters
    ==========
        - r,c : int
            - row and column-dimensions of the output.
        - k : int
            - diagonal offset.  In the output matrix,
            - out[m,m+k] == 1 for all m.
        - dtype : dtype
            - data-type of the output array.

    """
    warn("lil_eye is deprecated. use speye(... , format='lil') instead", \
            DeprecationWarning)
    return speye(r,c,k,dtype=dtype,format='lil')



#TODO remove this function
def lil_diags(diags,offsets,(m,n),dtype='d'):
    """Generate a lil_matrix with the given diagonals.

    Parameters
    ==========
        - diags : list of list of values e.g. [[1,2,3],[4,5]]
            - values to be placed on each indicated diagonal.
        - offsets : list of ints
            - diagonal offsets.  This indicates the diagonal on which
              the given values should be placed.
        - (r,c) : tuple of ints
            - row and column dimensions of the output.
        - dtype : dtype
            - output data-type.

    Example
    =======

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
    for k,diag in izip(offsets,diags):
        for ix,c in enumerate(xrange(clip(k,0,n),clip(m+k,0,n))):
            out.rows[c-k].append(c)
            out.data[c-k].append(diag[ix])
    return out


