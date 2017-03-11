"""
Routines for removing redundant (linearly dependent) equations from linear
programming equality constraints.
"""
# Author: Matt Haberland

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.linalg import svd


def _row_count(A):
    """
    Counts the number of nonzeros in each row of input array A.
    Nonzeros are defined as any element with absolute value greater than
    tol = 1e-13. This value should probably be an input to the function.

    Parameters
    ----------
    A : 2-D array
        An array representing a matrix

    Returns
    -------
    rowcount : 1-D array
        Number of nonzeros in each row of A

    """
    tol = 1e-13
    return np.sum(np.abs(A) > tol, axis=1)


def _get_densest(A, eligibleRows):
    """
    Returns the index of the densest row of A. Ignores rows that are not
    eligible for consideration.

    Parameters
    ----------
    A : 2-D array
        An array representing a matrix
    eligibleRows : 1-D logical array
        Values indicate whether the corresponding row of A is eligible
        to be considered

    Returns
    -------
    i_densest : int
        Index of the densest row in A eligible for consideration

    """
    rowCounts = _row_count(A)
    return np.argmax(rowCounts * eligibleRows)


def _remove_zero_rows(A, b):
    """
    Eliminates trivial equations from system of equations defined by Ax = b
   and identifies trivial infeasibilities

    Parameters
    ----------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations

    Returns
    -------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations
    status: int
        An integer indicating the status of the removal operation
        0: No infeasibility identified
        2: Trivially infeasible
    message : str
        A string descriptor of the exit status of the optimization.

    """
    status = 0
    message = ""
    i_zero = _row_count(A) == 0
    A = A[np.logical_not(i_zero), :]
    if not(np.allclose(b[i_zero], 0)):
        status = 2
        message = "There is a zero row in A_eq with a nonzero corresponding " \
                  "entry in b_eq. The problem is infeasible."
    b = b[np.logical_not(i_zero)]
    return A, b, status, message


def _remove_redundancy(A, b):
    """
    Eliminates redundant equations from system of equations defined by Ax = b
    and identifies infeasibilities.

    Parameters
    ----------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations

    Returns
    -------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations
    status: int
        An integer indicating the status of the system
        0: No infeasibility identified
        2: Trivially infeasible
    message : str
        A string descriptor of the exit status of the optimization.

    """

    A, b, status, message = _remove_zero_rows(A, b)

    if status != 0:
        return A, b, status, message

    U, s, Vh = svd(A)
    eps = np.finfo(float).eps
    tol = s.max() * max(A.shape) * eps

    m, n = A.shape
    s_min = s[-1] if m <= n else 0

    # this algorithm is inefficient
    # it relies on repeated singular value decomposition to find linearly
    # dependent rows (as identified by columns of U that correspond with zero
    # singular values). Unfortunately, only one row can be removed per
    # decomposition (I tried otherwise; doing so can cause problems.)
    # There are better algorithms out there such as:
    # Andersen, Erling D. "Finding all linearly dependent rows in
    # large-scale linear programming." Optimization Methods and Software
    # 6.3 (1995): 219-227.
    # but I was unable to get this to work.
    # It would be nice if we could do truncated SVD like sp.sparse.linalg.svds
    # but that function doesn't work well for the smallest singular value.

    while abs(s_min) < tol:
        v = U[:, -1]
        # rows need to be represented in significant amount
        eligibleRows = np.abs(v) > tol * 10e6
        if not np.any(eligibleRows) or np.any(np.abs(v.dot(A)) > tol):
            status = 4
            message = "Due to numerical issues, redundant equality " \
                      "constraints could not be removed automatically."
            break
        if np.any(np.abs(v.dot(b)) > tol):
            status = 2
            message = "There is a linear combination of rows of A_eq that " \
                      "results in zero, suggesting a redundant constraint. " \
                      "However the same linear combination of b_eq is " \
                      "nonzero, suggesting that the constraints conflict " \
                      "and the problem is infeasible."
            break

        i_remove = _get_densest(A, eligibleRows)
        A = np.delete(A, i_remove, axis=0)
        b = np.delete(b, i_remove)
        U, s, Vh = svd(A)
        m, n = A.shape
        s_min = s[-1] if m <= n else 0

    return A, b, status, message
