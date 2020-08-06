import numpy as np
from .misc import LinAlgError

__all__ = ['sinkhorn']


def sinkhorn(a):
    """
    Compute Sinkhorn normal form of a matrix.

    Calculate the Sinkhorm normal form ``A = D1 S D2`` where S is doubly
    stochastic and D1 and D2 are diagonal.

    Parameters
    ----------
    a : (N, N) array_like
        Matrix to be decomposed. Must be square and real with strictly positive
        elements.

    Returns
    -------
    D1 : ndarray
        Diagonal matrix of shape (N, N).
    S : ndarray
        Doubly stochastic matrix of shape (M, N).
    D2 : ndarray
        Diagonal matrix of shape (N, N).

    Raises
    ------
    LinAlgError
        Raised if normal form calculation fails

    Notes
    -----
    .. versionadded:: 1.6.0

    Examples
    --------
    >>> from scipy import linalg
    >>> a = np.random.uniform(size=(5, 5))

    >>> D1, S, D2 = linalg.sinkhorn(a)
    >>> np.allclose(a, D1 @ S @ D2)
    True

    """
    A = np.asarray_chkfinite(a)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("input matrix must be square")

    if np.any(A <= 0):
        raise ValueError("input matrix must be strictly positive")

    if np.any(np.iscomplex(A)):
        raise ValueError("input matrix must be real")

    N = len(A)
    r = np.ones(N)
    c = np.ones(N)
    S = A.copy()

    eps = np.finfo(a.dtype).eps * N
    for it in range(10000):
        prow = np.allclose(np.sum(S, axis=1), 1, atol=0, rtol=eps)
        pcol = np.allclose(np.sum(S, axis=0), 1, atol=0, rtol=eps)
        if prow and pcol:
            D1 = np.diag(1 / r)
            D2 = np.diag(1 / c)
            return D1, S, D2

        c = 1 / (r @ A)
        r = 1 / (A @ c)
        S = r[:, None] * A * c

    raise LinAlgError("Sinkhorn-Knopp did not converge.")
