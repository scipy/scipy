""" Sketching-based Matrix Computations """

# Author: Jordi Montes <jomsdev@gmail.com>
# August 28, 2017

from __future__ import division, print_function, absolute_import

import numpy as np


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

