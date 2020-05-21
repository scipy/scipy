import numpy as np
from . import linear_sum_assignment, OptimizeResult


def umeyama_eigendecomposition(A, B, maximize=False):
    """Finds an approximate solution to QAP using Umeyama's method."""
    # 'Correct' non-hermitian matrices
    if (A != A.T).any() or (B != B.T).any():
        A = (A + A.T) / 2 + 1j * (A - A.T) / 2
        B = (B + B.T) / 2 + 1j * (B - B.T) / 2

    va = np.linalg.eigh(A)[1]
    vb = np.linalg.eigh(B)[1]
    W = np.abs(vb) @ np.abs(va).T
    _, col_ind = linear_sum_assignment(W, maximize=maximize)

    # faster version of `float(np.trace(A @ B[col_ind][:, col_ind].T))`
    score = float(np.real(np.sum(A * B[col_ind][:, col_ind])))
    return OptimizeResult({"col_ind": col_ind, "score": score})


def quadratic_assignment(A, B, maximize=False):
    """Find an approximate solution to the quadratic assignment problem.

    This function finds approximate solutions to the quadratic assignment
    problem (QAP) [1]_, which is a problem of the form:

    .. math::

        \\min_P \\& \\ {\\text{trace}(APB^T P^T)}\\\\
        \\mbox{s.t.} \\& {P \\ \\epsilon \\ \\mathcal{P}}\\\\

    Here :math:`A` and :math:`B` are non-negative square matrices and
    :math:`\\mathcal{P}` is the set of all permutation matrices.

    The QAP is a NP-hard problem. This function uses Umeyama's
    eigendecomposition approach [2]_ to find an approximate solution. The
    results found are not guaranteed to be optimal.

    Parameters
    ----------
    A : 2d-array, square, non-negative
        A square matrix, also called the *cost* matrix.

    B : 2d-array, square, non-negative
        A square adjacency matrix, also called the *distance* matrix.

    maximize : bool (default: False)
        Maximizes the objective function if true.

    Returns
    -------
    res : OptimizeResult
        A :class:`scipy.optimize.OptimizeResult` consisting of the fields:
            col_ind : array
                An array of column indices giving the assignment.
            score : float
                The objective value of the assignment.

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Quadratic_assignment_problem

    .. [2] S Umeyama. An eigendecomposition approach to weighted graph matching
           problems. *IEEE Transactions on Pattern Analysis and Machine
           Intelligence* 10(5):695 - 703, September 1988,
           https://doi.org/10.1109/34.6778

    Examples
    --------

    >>> A = np.array([[0, 5, 8, 6],
    ...               [5, 0, 5, 1],
    ...               [8, 5, 0, 2],
    ...               [6, 1, 2, 0]])
    >>> B = np.array([[0, 1, 8, 4],
    ...               [1, 0, 5, 2],
    ...               [8, 5, 0, 5],
    ...               [4, 2, 5, 0]])
    >>> from scipy.optimize import quadratic_assignment
    >>> res = quadratic_assignment(A, B)
    >>> print(res)
     col_ind: array([3, 2, 1, 0])
       score: 200.0

    The `score` value, :math:`trace(APB^T P^T)`, can be calculated by forming
    the permutation matrix:
    >>> n = len(A)
    >>> P = np.zeros((n, n), dtype=int)
    >>> P[np.arange(n), res.col_ind] = 1
    >>> score = np.trace(A @ P @ B.T @ P.T)
    >>> print(score)
    200.0

    Alternatively, the score can be calculated by applying the appropriate
    permutation of ``B``:
    >>> score = np.trace(A @ B[res.col_ind][:, res.col_ind].T)
    >>> print(score)
    200.0
    """
    A = np.asarray(A)
    B = np.asarray(B)

    if A.ndim != 2:
        raise ValueError("``A`` must be a square matrix")
    if B.ndim != 2:
        raise ValueError("``B`` must be a square matrix")
    if A.shape != B.shape:
        raise ValueError("``A`` and ``B`` must be of equal size")
    if (A < 0).any():
        raise ValueError("``A`` contains negative entries")
    if (B < 0).any():
        raise ValueError("``B`` contains negative entries")

    return umeyama_eigendecomposition(A, B, maximize)
