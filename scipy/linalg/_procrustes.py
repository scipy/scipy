"""
Solve the orthogonal Procrustes problem.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.linalg import svd


__all__ = ['orthogonal_procrustes']


def orthogonal_procrustes(A, B, check_finite=True):
    """
    Compute the matrix solution of the orthogonal Procrustes problem.

    Given matrices A and B of equal shape, find an orthogonal matrix R
    that most closely maps A to B [1]_.
    Note that unlike higher level Procrustes analyses of spatial data,
    this function only uses orthogonal transformations like rotations
    and reflections, and it does not use scaling or translation.

    Parameters
    ----------
    A : (M, N) array_like
        Matrix to be mapped.
    B : (M, N) array_like
        Target matrix.

    Returns
    -------
    R : (N, N) ndarray
        The matrix solution of the orthogonal Procrustes problem.
        Minimizes the Frobenius norm of dot(A, R) - B, subject to
        dot(R.T, R) == I.

    References
    ----------
    .. [1] Peter H. Schonemann, "A generalized solution of the orthogonal
           Procrustes problem", Psychometrica -- Vol. 31, No. 1, March, 1996.

    """
    if check_finite:
        A = np.asarray_chkfinite(A)
        B = np.asarray_chkfinite(B)
    else:
        A = np.ascontiguousarray(A)
        B = np.asfortranarray(B)
    if A.ndim != 2:
        raise ValueError('expected ndim to be 2, but observed %s' % A.ndim)
    if A.shape != B.shape:
        raise ValueError('the shapes of A and B differ (%s vs %s)' % (
            A.shape, B.shape))
    u, w, vt = svd(A.T.dot(B))
    return u.dot(vt)

