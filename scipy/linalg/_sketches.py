""" Sketching-based Matrix Computations """

# Author: Jordi Montes <jomsdev@gmail.com>
# August 28, 2017

import numpy as np

from scipy._lib._util import (check_random_state, rng_integers,
                              _transition_to_rng, _apply_over_batch)

__all__ = ['clarkson_woodruff_transform']


def cwt_matrix(n_rows, n_columns, rng=None):
    r"""
    Generate a matrix S which represents a Clarkson-Woodruff transform.

    Given the desired size of matrix, the method returns a matrix S of size
    (n_rows, n_columns) where each column has all the entries set to 0
    except for one position which has been randomly set to +1 or -1 with
    equal probability.

    Parameters
    ----------
    n_rows : int
        Number of rows of S
    n_columns : int
        Number of columns of S
    rng : `numpy.random.Generator`, optional
        Pseudorandom number generator state. When `rng` is None, a new
        `numpy.random.Generator` is created using entropy from the
        operating system. Types other than `numpy.random.Generator` are
        passed to `numpy.random.default_rng` to instantiate a ``Generator``.

    Returns
    -------
    S : (n_rows, n_columns) csc_array
        The returned matrix has ``n_columns`` nonzero entries.

    Notes
    -----
    Given a matrix A, with probability at least 9/10,
    .. math:: \|SA\| = (1 \pm \epsilon)\|A\|
    Where the error epsilon is related to the size of S.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import cwt_matrix
    >>> rng = np.random.default_rng(seed=42)
    >>> S = cwt_matrix(3, 6, rng=rng)
    >>> S.toarray()
    array([[ 1,  0,  0,  0,  0,  0],
           [ 0,  0,  1,  1, -1,  0],
           [ 0, -1,  0,  0,  0, -1]])

    The matrix has exactly one nonzero entry per column, randomly set to
    +1 or -1:

    >>> S.nnz == 6  # one nonzero per column
    True

    """
    # lazy import to prevent to prevent sparse dependency for whole module (gh-23420)
    from scipy.sparse import csc_array
    rng = check_random_state(rng)
    rows = rng_integers(rng, 0, n_rows, n_columns)
    cols = np.arange(n_columns+1)
    signs = rng.choice([1, -1], n_columns)
    S = csc_array((signs, rows, cols), shape=(n_rows, n_columns))
    return S