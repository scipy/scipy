"""Functions to construct sparse matrices
"""

__docformat__ = "restructuredtext en"

__all__ = [ 'spdiags', 'eye', 'identity', 'kron', 'kronsum',
            'hstack', 'vstack', 'bmat', 'rand']


from warnings import warn

import numpy as np

from sputils import upcast

from csr import csr_matrix
from csc import csc_matrix
from bsr import bsr_matrix
from coo import coo_matrix
from lil import lil_matrix
from dia import dia_matrix

def spdiags(data, diags, m, n, format=None):
    """
    Return a sparse matrix from diagonals.

    Parameters
    ----------
    data   : array_like
        matrix diagonals stored row-wise
    diags  : diagonals to set
        - k = 0  the main diagonal
        - k > 0  the k-th upper diagonal
        - k < 0  the k-th lower diagonal
    m, n : int
        shape of the result
    format : format of the result (e.g. "csr")
        By default (format=None) an appropriate sparse matrix
        format is returned.  This choice is subject to change.

    See Also
    --------
    dia_matrix : the sparse DIAgonal format.

    Examples
    --------
    >>> data = array([[1,2,3,4],[1,2,3,4],[1,2,3,4]])
    >>> diags = array([0,-1,2])
    >>> spdiags(data, diags, 4, 4).todense()
    matrix([[1, 0, 3, 0],
            [1, 2, 0, 4],
            [0, 2, 3, 0],
            [0, 0, 3, 4]])

    """
    return dia_matrix((data, diags), shape=(m,n)).asformat(format)

def identity(n, dtype='d', format=None):
    """Identity matrix in sparse format

    Returns an identity matrix with shape (n,n) using a given
    sparse format and dtype.

    Parameters
    ----------
    n : integer
        Shape of the identity matrix.
    dtype :
        Data type of the matrix
    format : string
        Sparse format of the result, e.g. format="csr", etc.

    Examples
    --------
    >>> identity(3).todense()
    matrix([[ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.]])
    >>> identity(3, dtype='int8', format='dia')
    <3x3 sparse matrix of type '<type 'numpy.int8'>'
            with 3 stored elements (1 diagonals) in DIAgonal format>

    """

    if format in ['csr','csc']:
        indptr  = np.arange(n+1, dtype=np.intc)
        indices = np.arange(n,   dtype=np.intc)
        data    = np.ones(n,     dtype=dtype)
        cls = eval('%s_matrix' % format)
        return cls((data,indices,indptr),(n,n))
    elif format == 'coo':
        row  = np.arange(n, dtype=np.intc)
        col  = np.arange(n, dtype=np.intc)
        data = np.ones(n, dtype=dtype)
        return coo_matrix((data,(row,col)),(n,n))
    elif format == 'dia':
        data = np.ones(n, dtype=dtype)
        diags = [0]
        return dia_matrix((data,diags), shape=(n,n))
    else:
        return identity(n, dtype=dtype, format='csr').asformat(format)


def eye(m, n, k=0, dtype='d', format=None):
    """eye(m, n) returns a sparse (m x n) matrix where the k-th diagonal
    is all ones and everything else is zeros.
    """
    m,n = int(m),int(n)
    diags = np.ones((1, max(0, min(m + k, n))), dtype=dtype)
    return spdiags(diags, k, m, n).asformat(format)


def kron(A, B, format=None):
    """kronecker product of sparse matrices A and B

    Parameters
    ----------
    A : sparse or dense matrix
        first matrix of the product
    B : sparse or dense matrix
        second matrix of the product
    format : string
        format of the result (e.g. "csr")

    Returns
    -------
    kronecker product in a sparse matrix format


    Examples
    --------
    >>> A = csr_matrix(array([[0,2],[5,0]]))
    >>> B = csr_matrix(array([[1,2],[3,4]]))
    >>> kron(A,B).todense()
    matrix([[ 0,  0,  2,  4],
            [ 0,  0,  6,  8],
            [ 5, 10,  0,  0],
            [15, 20,  0,  0]])

    >>> kron(A,[[1,2],[3,4]]).todense()
    matrix([[ 0,  0,  2,  4],
            [ 0,  0,  6,  8],
            [ 5, 10,  0,  0],
            [15, 20,  0,  0]])

    """
    B = coo_matrix(B)

    if (format is None or format == "bsr") and 2*B.nnz >= B.shape[0] * B.shape[1]:
        #B is fairly dense, use BSR
        A = csr_matrix(A,copy=True)

        output_shape = (A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])

        if A.nnz == 0 or B.nnz == 0:
            # kronecker product is the zero matrix
            return coo_matrix( output_shape )

        B = B.toarray()
        data = A.data.repeat(B.size).reshape(-1,B.shape[0],B.shape[1])
        data = data * B

        return bsr_matrix((data,A.indices,A.indptr), shape=output_shape)
    else:
        #use COO
        A = coo_matrix(A)
        output_shape = (A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])

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

def kronsum(A, B, format=None):
    """kronecker sum of sparse matrices A and B

    Kronecker sum of two sparse matrices is a sum of two Kronecker
    products kron(I_n,A) + kron(B,I_m) where A has shape (m,m)
    and B has shape (n,n) and I_m and I_n are identity matrices
    of shape (m,m) and (n,n) respectively.

    Parameters
    ----------
    A
        square matrix
    B
        square matrix
    format : string
        format of the result (e.g. "csr")

    Returns
    -------
    kronecker sum in a sparse matrix format

    Examples
    --------


    """
    A = coo_matrix(A)
    B = coo_matrix(B)

    if A.shape[0] != A.shape[1]:
        raise ValueError('A is not square')

    if B.shape[0] != B.shape[1]:
        raise ValueError('B is not square')

    dtype = upcast(A.dtype, B.dtype)

    L = kron(identity(B.shape[0],dtype=dtype), A, format=format)
    R = kron(B, identity(A.shape[0],dtype=dtype), format=format)

    return (L+R).asformat(format) #since L + R is not always same format


def hstack(blocks, format=None, dtype=None):
    """
    Stack sparse matrices horizontally (column wise)

    Parameters
    ----------
    blocks
        sequence of sparse matrices with compatible shapes
    format : string
        sparse format of the result (e.g. "csr")
        by default an appropriate sparse matrix format is returned.
        This choice is subject to change.

    See Also
    --------
    vstack : stack sparse matrices vertically (row wise)

    Examples
    --------
    >>> from scipy.sparse import coo_matrix, vstack
    >>> A = coo_matrix([[1,2],[3,4]])
    >>> B = coo_matrix([[5],[6]])
    >>> hstack( [A,B] ).todense()
    matrix([[1, 2, 5],
            [3, 4, 6]])

    """
    return bmat([blocks], format=format, dtype=dtype)

def vstack(blocks, format=None, dtype=None):
    """
    Stack sparse matrices vertically (row wise)

    Parameters
    ----------
    blocks
        sequence of sparse matrices with compatible shapes
    format : string
        sparse format of the result (e.g. "csr")
        by default an appropriate sparse matrix format is returned.
        This choice is subject to change.

    See Also
    --------
    hstack : stack sparse matrices horizontally (column wise)

    Examples
    --------
    >>> from scipy.sparse import coo_matrix, vstack
    >>> A = coo_matrix([[1,2],[3,4]])
    >>> B = coo_matrix([[5,6]])
    >>> vstack( [A,B] ).todense()
    matrix([[1, 2],
            [3, 4],
            [5, 6]])

    """
    return bmat([ [b] for b in blocks ], format=format, dtype=dtype)

def bmat(blocks, format=None, dtype=None):
    """
    Build a sparse matrix from sparse sub-blocks

    Parameters
    ----------
    blocks
        grid of sparse matrices with compatible shapes
        an entry of None implies an all-zero matrix
    format : sparse format of the result (e.g. "csr")
        by default an appropriate sparse matrix format is returned.
        This choice is subject to change.

    Examples
    --------
    >>> from scipy.sparse import coo_matrix, bmat
    >>> A = coo_matrix([[1,2],[3,4]])
    >>> B = coo_matrix([[5],[6]])
    >>> C = coo_matrix([[7]])
    >>> bmat( [[A,B],[None,C]] ).todense()
    matrix([[1, 2, 5],
            [3, 4, 6],
            [0, 0, 7]])

    >>> bmat( [[A,None],[None,C]] ).todense()
    matrix([[1, 2, 0],
            [3, 4, 0],
            [0, 0, 7]])

    """

    blocks = np.asarray(blocks, dtype='object')

    if np.rank(blocks) != 2:
        raise ValueError('blocks must have rank 2')

    M,N = blocks.shape

    block_mask   = np.zeros(blocks.shape,    dtype=np.bool)
    brow_lengths = np.zeros(blocks.shape[0], dtype=np.intc)
    bcol_lengths = np.zeros(blocks.shape[1], dtype=np.intc)

    # convert everything to COO format
    for i in range(M):
        for j in range(N):
            if blocks[i,j] is not None:
                A = coo_matrix(blocks[i,j])
                blocks[i,j] = A
                block_mask[i,j] = True

                if brow_lengths[i] == 0:
                    brow_lengths[i] = A.shape[0]
                else:
                    if brow_lengths[i] != A.shape[0]:
                        raise ValueError('blocks[%d,:] has incompatible row dimensions' % i)

                if bcol_lengths[j] == 0:
                    bcol_lengths[j] = A.shape[1]
                else:
                    if bcol_lengths[j] != A.shape[1]:
                        raise ValueError('blocks[:,%d] has incompatible column dimensions' % j)


    # ensure that at least one value in each row and col is not None
    if brow_lengths.min() == 0:
        raise ValueError('blocks[%d,:] is all None' % brow_lengths.argmin() )
    if bcol_lengths.min() == 0:
        raise ValueError('blocks[:,%d] is all None' % bcol_lengths.argmin() )

    nnz = sum([ A.nnz for A in blocks[block_mask] ])
    if dtype is None:
        dtype = upcast( *tuple([A.dtype for A in blocks[block_mask]]) )

    row_offsets = np.concatenate(([0], np.cumsum(brow_lengths)))
    col_offsets = np.concatenate(([0], np.cumsum(bcol_lengths)))

    data = np.empty(nnz, dtype=dtype)
    row  = np.empty(nnz, dtype=np.intc)
    col  = np.empty(nnz, dtype=np.intc)

    nnz = 0
    for i in range(M):
        for j in range(N):
            if blocks[i,j] is not None:
                A = blocks[i,j]
                data[nnz:nnz + A.nnz] = A.data
                row[nnz:nnz + A.nnz]  = A.row
                col[nnz:nnz + A.nnz]  = A.col

                row[nnz:nnz + A.nnz] += row_offsets[i]
                col[nnz:nnz + A.nnz] += col_offsets[j]

                nnz += A.nnz

    shape = (np.sum(brow_lengths), np.sum(bcol_lengths))
    return coo_matrix((data, (row, col)), shape=shape).asformat(format)

def rand(m, n, density=0.01, format="coo", dtype=None):
    """Generate a sparse matrix of the given shape and density with uniformely
    distributed values.

    Parameters
    ----------
    m, n: int
        shape of the matrix
    density: real
        density of the generated matrix: density equal to one means a full
        matrix, density of 0 means a matrix with no non-zero items.
    format: str
        sparse matrix format.
    dtype: dtype
        type of the returned matrix values.

    Notes
    -----
    Only float types are supported for now.
    """
    if density < 0 or density > 1:
        raise ValueError("density expected to be 0 <= density <= 1")
    if dtype and not dtype in [np.float32, np.float64, np.longdouble]:
        raise NotImplementedError("type %s not supported" % dtype)

    mn = m * n

    # XXX: sparse uses intc instead of intp...
    tp = np.intp
    if mn > np.iinfo(tp).max:
        msg = """\
Trying to generate a random sparse matrix such as the product of dimensions is
greater than %d - this is not supported on this machine
"""
        raise ValueError(msg % np.iinfo(tp).max)

    # Number of non zero values
    k = long(density * m * n)

    # Generate a few more values than k so that we can get unique values
    # afterwards.
    # XXX: one could be smarter here
    mlow = 5
    fac = 1.02
    gk = min(k + mlow, fac * k)

    def _gen_unique_rand(_gk):
        id = np.random.rand(_gk)
        return np.unique(np.floor(id * mn))[:k]

    id = _gen_unique_rand(gk)
    while id.size < k:
        gk *= 1.05
        id = _gen_unique_rand(gk)

    j = np.floor(id * 1. / m).astype(tp)
    i = (id - j * m).astype(tp)
    vals = np.random.rand(k).astype(dtype)
    return coo_matrix((vals, (i, j)), shape=(m, n)).asformat(format)
