""" Sketching-based Matrix Computations """

# Author: Jordi Montes <jomsdev@gmail.com>
# August 28, 2017

from __future__ import division, print_function, absolute_import

import numpy as np

__all__ = ['clarkson_woodruff_transform']


def cwt_matrix(n_rows, n_columns):
    """
    Generate a matrix S to be used in the Clarkson-Woodruff sketch.
    Given the size of a desired matrix, the method returns a matrix S of size (n_rows, n_columns)
    where each column has all the entries set to 0 less one position with has been randomly set to
    +1 o -1 with equal probability.

    Parameters
    ----------
    n_rows: int
        number of rows of S
    n_columns: int
        number of columns of S

    Returns
    -------
    S : (n_rows, n_columns) array_like

    Notes
    -----
    Given a matrix A, with probability at least 9/10, norm(SA) == (1+epsilon)*norm(A)
    Where epsilon is related to the size of S
    """
    S = np.zeros((n_rows, n_columns))
    nz_positions = np.random.randint(0, n_rows, n_columns)
    values = np.random.choice([1, -1], n_columns)
    for i in range(n_columns):
        S[nz_positions[i]][i] = values[i]

    return S

def clarkson_woodruff_transform(input_matrix, sketch_size):
    """
    Given an input_matrix  of size (n, d), compute a matrix A' of size  (sketch_size, d)
    which holds:
    $||Ax|| = (1 \pm \epsilon) ||A'x||$
    with high probability.

    The error is related to the number of rows of the sketch. sketch_size = $poly(r*\epsilon^{-1})$
    
    Parameters
    ----------
    input_matrix: (n, d) array_like
        Input matrix
    sketch_size: int
        number of rows for the sketch

    Returns
    -------
    A' : (sketch_size, d) array_like
        Sketch of A

    Notes
    -----
    This is an implementation of the Clarckson-Woodruff Transform (also known as CountSketch) introduced for
    first time in Kenneth L. Clarkson and David P. Woodruff. Low rank approximation and regression in input sparsity time. In STOC, 2013.
    A' can be computed in O(nnz(A)) but we don't take advantage of sparse matrix in this implementation
    """

    S = cwt_matrix(sketch_size, input_matrix.shape[0])
    return np.dot(S, input_matrix)
