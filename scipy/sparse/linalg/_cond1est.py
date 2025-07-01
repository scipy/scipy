from scipy.sparse.linalg import norm, splu

__all__ = ['invnormest', 'cond1est']


def invnormest(A, ord=None):
    """Compute an estimate for the norm of the inverse of a sparse matrix.

    Parameters
    ----------
    A : sparray or LinearOperator
        A square, sparse matrix. Any matrix not in CSC format will be converted
        internally, and raise a SparseEfficiencyWarning.
    ord : {1, inf}, optional
        Order of the norm. If None, defaults to the 1-norm. inf means numpy's
        `inf` object.

    Returns
    -------
    est : float
        An estimate of the norm of the matrix inverse.

    Notes
    -----
    This function computes an LU decomposition using SuperLU, and then runs the
    appropriate `gscon` procedure for the data type. Use `SuperLU.invnormest`
    if you already have an LU decomposition.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import invnormest
    >>> A = csc_matrix([[1., 0., 0.], [5., 8., 2.], [0., -1., 0.]], dtype=float)
    >>> A.toarray()
    array([[ 1.,  0.,  0.],
           [ 5.,  8.,  2.],
           [ 0., -1.,  0.]])
    >>> invnormest(A, ord=1)
    5.0
    >>> float(np.linalg.norm(np.linalg.inv(A.toarray()), ord=1))
    5.0
    """
    lu = splu(A)
    return lu.invnormest(ord=ord)


def cond1est(A):
    """
    Compute an estimate for the reciprocal of the condition number 
    of a sparse matrix, using 1-norms.

    Parameters
    ----------
    A : ndarray or other linear operator
        A sparse matrix for which an LU matrix can be computed. CSC would
        be most efficient.

    Returns
    -------
    cond : float
        An estimate of the 1-norm condition number of A.

    Notes
    -----
    Computes an LU decomposition and runs the gscon procedure from SuperLU.
    Use scipy.sparse.linalg.SuperLU.invnormest if you already have an LU 
    decomposition.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import cond1est
    >>> A = csc_matrix([[1., 0., 0.], [5., 8., 2.], [0., -1., 0.]], dtype=float)
    >>> A.toarray()
    array([[ 1.,  0.,  0.],
           [ 5.,  8.,  2.],
           [ 0., -1.,  0.]])
    >>> cond1est(A)
    45.0
    >>> np.linalg.cond(A.toarray(), p=1)
    45.0
    """
    return norm(A, 1) * invnormest(A, ord=1)
