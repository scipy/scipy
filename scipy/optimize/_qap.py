import numpy as np
from . import linear_sum_assignment


def quadratic_assignment(A, B, maximize=False):
    """Find an approximate solution to the quadratic assignment problem.

    This function finds approximate solutions to the quadratic assignment
    problem (QAP), which is a problem of the form:

    .. math::

        \\min_P \\& \\ {-\\text{trace}(APB^T P^T)}\\\\
        \\mbox{s.t. } \\& {P \\ \\epsilon \\ \\mathcal{P}}\\\\

    Here and :math:`A` and :math:`B` are non-negative square matrices and
    :math:`\\mathcal{P}` is the set of all permutation matrices.

    The QAP is a NP-hard problem. This function uses Umeyama's
    eigendecomposition approach to find an approximate solution. The results
    found are not guaranteed to be optimal.

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
    row_ind, col_ind : array
        An array of row indices and one of corresponding column indices giving
        the assignment. The cost of the assignment can be computed as
        ``np.trace(A.T[col_ind] @ B[col_ind].T)``. The row indices will be
        sorted.

    References
    ----------

    1. https://en.wikipedia.org/wiki/Quadratic_assignment_problem

    2. S Umeyama. An eigendecomposition approach to weighted graph matching problems.
       *IEEE Transactions on Pattern Analysis and Machine Intelligence*
       10(5):695 - 703, September 1988, https://doi.org/10.1109/34.6778

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
    >>> row_ind, col_ind = quadratic_assignment(A, B)
    >>> print(row_ind)
    array([0, 1, 2, 3])
    >>> print(col_ind)
    array([3, 2, 1, 0])
    """
    A = np.asarray(A)
    B = np.asarray(B)

    if A.ndim != 2:
        raise ValueError("``A`` must be a square matrix")
    if B.ndim != 2:
        raise ValueError("``B`` must be a square matrix")
    if A.shape != B.shape:
        raise ValueError("A and B must be of equal size")
    if (A < 0).any():
        raise ValueError("``A`` contains negative entries")
    if (B < 0).any():
        raise ValueError("``B`` contains negative entries")

    def eigenvectors(x):
        # get eigenvectors sorted by eigenvalue magnitude (large to small)
        l, v = np.linalg.eig(x)
        indices = np.argsort(l, kind='merge')[::-1]
        return v[:, indices]

    cost_matrix = np.abs(eigenvectors(B)) @ np.abs(eigenvectors(A)).T
    if maximize:
        cost_matrix = -cost_matrix
    return linear_sum_assignment(cost_matrix)
