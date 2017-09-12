""" Sketching-based Matrix Computations """

# Author: Jordi Montes <jomsdev@gmail.com>
# August 28, 2017

from __future__ import division, print_function, absolute_import

import numpy as np

from scipy._lib._util import check_random_state

__all__ = ['clarkson_woodruff_transform']


def cwt_matrix(n_rows, n_columns, seed=None):
    r""""
    Generate a matrix S for the Clarkson-Woodruff sketch. Given the desired
    size of matrix, the method returns a matrix S of size (n_rows, n_columns)
    where each column has all the entries set to 0 less one position which has
    been randomly set to +1 o -1 with equal probability.

    Parameters
    ----------
    n_rows: int
        number of rows of S
    n_columns: int
        number of columns of S
    seed: int
        Seed the generator

    Returns
    -------
    S : (n_rows, n_columns) array_like

    Notes
    -----
    Given a matrix A, with probability at least 9/10,
    .. math:: ||SA|| == (1 \pm \epsilon)||A||
    Where epsilon is related to the size of S
    """
    S = np.zeros((n_rows, n_columns))
    nz_positions = np.random.randint(0, n_rows, n_columns)
    rng = check_random_state(seed)
    values = rng.choice([1, -1], n_columns)
    for i in range(n_columns):
        S[nz_positions[i]][i] = values[i]

    return S

def clarkson_woodruff_transform(input_matrix, sketch_size, seed=None):
    r""""
    Given an input_matrix  of size (n, d), compute a matrix A' of size
    (sketch_size, d) which holds:
    
    .. math:: ||Ax|| = (1 \pm \epsilon)||A'x||

    with high probability.

    The error is related to the number of rows of the sketch and it is bounded
    
    .. math:: poly(r(\epsilon^{-1}))
    
    Parameters
    ----------
    input_matrix: (n, d) array_like
        Input matrix
    sketch_size: int
        number of rows for the sketch
    seed: int
        Seed the generator

    Returns
    -------
    A' : (sketch_size, d) array_like
        Sketch of A

    Notes
    -----
    This is an implementation of the Clarkson-Woodruff Transform (CountSketch)
    A' can be computed in O(nnz(A)) but we don't take advantage of sparse
    matrices in this implementation

    Examples
    --------
    Given a big dense matrix `A`:

    >>> from scipy import linalg
    >>> n_rows, n_columns, sketch_n_rows = (2000, 100, 100)
    >>> threshold = 0.1
    >>> tmp = np.random.normal(0, 0.1, n_rows*n_columns)
    >>> A = np.reshape(tmp, (n_rows, n_columns))
    >>> sketch = linalg.clarkson_woodruff_transform(A, sketch_n_rows)
    >>> sketch.shape
    (100, 100)
    >>> normA = linalg.norm(A)
    >>> normSketch = linalg.norm(sketch)
    >>> # with high probability
    >>> # abs(normA-normSketch) < threshold

    References
    ----------
    .. [1] Kenneth L. Clarkson and David P. Woodruff. Low rank approximation and
           regression in input sparsity time. In STOC, 2013.
    """

    S = cwt_matrix(sketch_size, input_matrix.shape[0], seed)
    return np.dot(S, input_matrix)
