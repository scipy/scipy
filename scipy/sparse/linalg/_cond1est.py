import numpy as np

from scipy.sparse import issparse
from scipy.sparse.linalg import splu

__all__ = ['cond1est']


def cond1est(A):
    r"""Estimate the condition number of a sparse matrix in the 1-norm.

    This function computes an LU decomposition of ``A`` and calls
    ``SuperLU.cond1est``. If you already have an LU factorization computed, use
    ``SuperLU.cond1est`` directly.

    Parameters
    ----------
    A : (N, N) sparray
        A square, sparse matrix. Any matrix not in CSC format will be converted
        internally, and raise a ``SparseEfficiencyWarning``.

    Returns
    -------
    cond : float
        An estimate of the condition number of ``A`` in the 1-norm.

    See Also
    --------
    numpy.linalg.cond : Compute the condition number of a dense matrix.
    SuperLU.cond1est : Compute an estimate of the condition number of a sparse
                       matrix.

    Notes
    -----
    The condition number is defined as [0]_:

    .. math:: \kappa(A) = \left\| A \right\|_1 \left\| A^{-1} \right\|_1.

    This function computes the 1-norm of the matrix exactly, and a lower bound
    estimate for the 1-norm of the inverse using SuperLU internals. It is
    similar to ``np.linalg.cond(A, p=1)`` for dense matrices, but given that
    this function uses an estimate, results will not be identical, especially
    for ill-conditioned matrices.

    .. versionadded:: 1.17.0

    References
    ----------
    .. [0] "Condition Number", Wikipedia,
           https://en.wikipedia.org/wiki/Condition_number

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import cond1est
    >>> A = csc_array([[1., 0., 0.], [5., 8., 2.], [0., -1., 0.]], dtype=float)
    >>> A.toarray()
    array([[ 1.,  0.,  0.],
           [ 5.,  8.,  2.],
           [ 0., -1.,  0.]])
    >>> cond1est(A)
    np.float64(45.0)
    >>> np.linalg.cond(A.toarray(), p=1)
    np.float64(45.0)
    """
    if not issparse(A):
       raise TypeError("Input is not a sparse matrix. "
                       "Use np.linalg.cond for dense matrices.")

    if A.ndim != 2:
        raise ValueError("Input must be a 2-dimensional matrix.")

    M, N = A.shape

    if M != N:
        raise ValueError("Matrix must be square.")

    if M == 0:
        raise ValueError("Condition number of an empty matrix is undefined.")

    try:
        return splu(A).cond1est()
    except RuntimeError as e:
        if "Factor is exactly singular" in str(e):
            return np.inf
        else:
            raise e
