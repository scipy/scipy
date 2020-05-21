import numpy as np
from . import linear_sum_assignment, OptimizeResult


def eigenvectors(x):
    # get eigenvectors sorted by eigenvalue magnitude (large to small)
    l, v = np.linalg.eig(x)
    indices = np.argsort(l, kind='merge')[::-1]
    return v[:, indices]


def umeyama_eigendecomposition(A, B, maximize=False):
    # Finds an approximate solution to QAP using Umeyama's method
    W = np.abs(eigenvectors(B)) @ np.abs(eigenvectors(A)).T
    _, col_ind = linear_sum_assignment(W, maximize=maximize)
    score = float(np.trace(A @ B[col_ind][:, col_ind].T))
    return OptimizeResult({"col_ind": col_ind, "score": score})
