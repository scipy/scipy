import numpy as np

from scipy.sparse.linalg import norm, splu

__all__ = ['normest_inv', 'cond1est']


def normest_inv(A, ord=None):
    """Compute an estimate of the norm of the inverse of a sparse matrix.

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

    See Also
    --------
    norm : Compute the norm of a sparse matrix.
    numpy.linalg.norm : Compute the norm of a dense matrix.
    cond1est : Compute an estimate for the reciprocal of the condition number
        of a sparse matrix.
    splu : Compute the LU decomposition of a sparse matrix.

    Notes
    -----
    This function computes an LU decomposition using SuperLU, and then runs the
    appropriate ``gscon``[0]_ procedure for the data type. Use
    ``SuperLU.normest_inv`` if you already have an LU decomposition.

    .. versionadded:: 1.17.0

    References
    ----------
    .. [0] https://portal.nersc.gov/project/sparse/superlu/superlu_code_html/dgscon_8c.html

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import normest_inv
    >>> A = csc_array([[1., 0., 0.], [5., 8., 2.], [0., -1., 0.]], dtype=float)
    >>> A.toarray()
    array([[ 1.,  0.,  0.],
           [ 5.,  8.,  2.],
           [ 0., -1.,  0.]])
    >>> normest_inv(A, ord=1)
    5.0
    >>> np.linalg.norm(np.linalg.inv(A.toarray()), ord=1)
    np.float64(5.0)
    """
    return splu(A).normest_inv(ord=ord)


def cond1est(A):
    r"""Estimate the condition number of a sparse matrix in the 1-norm.

    Parameters
    ----------
    A : ndarray or other linear operator
        A square, sparse matrix. Any matrix not in CSC format will be converted
        internally, and raise a SparseEfficiencyWarning.

    Returns
    -------
    cond : float
        An estimate of the condition number of A in the 1-norm.

    See Also
    --------
    normest_inv : Compute an estimate for the norm of the matrix inverse.
    numpy.linalg.cond : Compute the condition number of a dense matrix.

    Notes
    -----
    The condition number is defined as[0]_:

    .. math:: \kappa(A) = \left\| A \right\|_1 \left\| A^{-1} \right\|_1.

    This function computes the 1-norm of the matrix and a lower bound estimate
    for the 1-norm of the inverse using ``normest_inv``. It is similar to 
    ``np.linalg.cond(A, p=1)`` for dense matrices, but given that this function
    uses an estimate, results will not be identical, especially for
    ill-conditioned matrices.

    .. versionadded:: 1.17.0

    References
    ----------
    .. [0] https://en.wikipedia.org/wiki/Condition_number

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
    np.float64(45.0)
    >>> np.linalg.cond(A.toarray(), p=1)
    np.float64(45.0)
    """
    M, N = A.shape

    if M != N:
        raise ValueError("Matrix must be square.")

    if M == 0:
        raise ValueError("Condition number of an empty matrix is undefined.")

    # Compute the 1-norm of A exactly
    norm_A = norm(A, 1)

    try:
        norm_A_inv = normest_inv(A, ord=1)
    except RuntimeError as e:
        if "Factor is exactly singular" in str(e):
            return np.inf
        else:
            raise e

    return norm_A * norm_A_inv
