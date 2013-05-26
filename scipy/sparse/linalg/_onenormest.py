"""Sparse block 1-norm estimator.
"""

from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
import scipy.sparse.linalg
from scipy.sparse.linalg import aslinearoperator


__all__ = ['onenormest']


def onenormest(A, t=2, itmax=5, compute_v=False, compute_w=False):
    """
    Compute a lower bound of the 1-norm of a sparse matrix.

    .. versionadded:: 0.13.0

    Parameters
    ----------
    A : ndarray or other linear operator
        A linear operator that can be transposed and that can
        produce matrix products.
    t : int, optional
        A positive parameter controlling the tradeoff between
        accuracy versus time and memory usage.
        Larger values take longer and use more memory
        but give more accurate output.
    itmax : int, optional
        Use at most this many iterations.
    compute_v : bool, optional
        Request a norm-maximizing linear operator input vector if True.
    compute_w : bool, optional
        Request a norm-maximizing linear operator output vector if True.

    Returns
    -------
    est : float
        An underestimate of the 1-norm of the sparse matrix.
    v : ndarray, optional
        The vector such that ||Av||_1 == est*||v||_1.
        It can be thought of as an input to the linear operator
        that gives an output with particularly large norm.
    w : ndarray, optional
        The vector Av which has relatively large 1-norm.
        It can be thought of as an output of the linear operator
        that is relatively large in norm compared to the input.

    Notes
    -----
    This is algorithm 2.4 of [1].

    In [2] it is described as follows.
    "This algorithm typically requires the evaluation of
    about 4t matrix-vector products and almost invariably
    produces a norm estimate (which is, in fact, a lower
    bound on the norm) correct to within a factor 3."

    References
    ----------
    .. [1] Nicholas J. Higham and Francoise Tisseur (2000),
           "A Block Algorithm for Matrix 1-Norm Estimation,
           with an Application to 1-Norm Pseudospectra."
           SIAM J. Matrix Anal. Appl. Vol. 21, No. 4, pp. 1185-1201.

    .. [2] Awad H. Al-Mohy and Nicholas J. Higham (2009),
           "A new scaling and squaring algorithm for the matrix exponential."
           SIAM J. Matrix Anal. Appl. Vol. 31, No. 3, pp. 970-989.

    """

    # Check the input.
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected the operator to act like a square matrix')

    # If the operator size is small compared to t,
    # then it is easier to compute the exact norm.
    # Otherwise estimate the norm.
    n = A.shape[1]
    if t >= n:
        A_explicit = np.asarray(aslinearoperator(A).matmat(np.identity(n)))
        if A_explicit.shape != (n, n):
            raise Exception('internal error: ',
                    'unexpected shape ' + str(A_explicit.shape))
        col_abs_sums = abs(A_explicit).sum(axis=0)
        if col_abs_sums.shape != (n, ):
            raise Exception('internal error: ',
                    'unexpected shape ' + str(col_abs_sums.shape))
        argmax_j = np.argmax(col_abs_sums)
        v = elementary_vector(n, argmax_j)
        w = A_explicit[:, argmax_j]
        est = col_abs_sums[argmax_j]
    else:
        est, v, w, nmults, nresamples = _onenormest_core(A, A.T, t, itmax)

    # Report the norm estimate along with some certificates of the estimate.
    if compute_v or compute_w:
        result = (est,)
        if compute_v:
            result += (v,)
        if compute_w:
            result += (w,)
        return result
    else:
        return est


def sign_round_up(X):
    # differs from numpy sign in the way it handles entries that are zero
    return np.sign(np.sign(X) + 0.5)


def elementary_vector(n, i):
    v = np.zeros(n, dtype=float)
    v[i] = 1
    return v


def vectors_are_parallel(v, w):
    # Columns are considered parallel when they are equal or negative.
    # Entries are required to be in {-1, 1},
    # which guarantees that the magnitudes of the vectors are identical.
    if v.ndim != 1 or v.shape != w.shape:
        raise ValueError('expected conformant vectors with entries in {-1,1}')
    n = v.shape[0]
    return np.dot(v, w) == n


def every_col_of_X_is_parallel_to_a_col_of_Y(X, Y):
    for v in X.T:
        if not any(vectors_are_parallel(v, w) for w in Y.T):
            return False
    return True


def column_needs_resampling(i, X, Y=None):
    # column i of X needs resampling if either
    # it is parallel to a previous column of X or
    # it is parallel to a column of Y
    n, t = X.shape
    v = X[:, i]
    if any(vectors_are_parallel(v, X[:, j]) for j in range(i)):
        return True
    if Y is not None:
        if any(vectors_are_parallel(v, w) for w in Y.T):
            return True
    return False


def resample_column(i, X):
    X[:, i] = np.random.randint(0, 2, size=X.shape[0])*2 - 1


def norm_1d_1(v):
    # this is faster than calling the numpy/scipy norm function
    return np.sum(np.abs(v))


def norm_1d_inf(v):
    # this is faster than calling the numpy/scipy norm function
    return np.max(np.abs(v))


def less_than_or_close(a, b):
    return np.allclose(a, b) or (a < b)


def _algorithm_2_2(A, AT, t):
    """
    This is Algorithm 2.2.

    Parameters
    ----------
    A : ndarray or other linear operator
        A linear operator that can produce matrix products.
    AT : ndarray or other linear operator
        The transpose of A.
    t : int, optional
        A positive parameter controlling the tradeoff between
        accuracy versus time and memory usage.

    Returns
    -------
    g : sequence
        A non-negative decreasing vector
        such that g[j] is a lower bound for the 1-norm
        of the column of A of jth largest 1-norm.
        The first entry of this vector is therefore a lower bound
        on the 1-norm of the linear operator A.
        This sequence has length t.
    ind : sequence
        The ith entry of ind is the index of the column A whose 1-norm
        is given by g[i].
        This sequence of indices has length t, and its entries are
        chosen from range(n), possibly with repetition,
        where n is the order of the operator A.

    Notes
    -----
    This algorithm is mainly for testing.
    It uses the 'ind' array in a way that is similar to
    its usage in algorithm 2.4.  This algorithm 2.2 may be easier to test,
    so it gives a chance of uncovering bugs related to indexing
    which could have propagated less noticably to algorithm 2.4.

    """
    A_linear_operator = aslinearoperator(A)
    AT_linear_operator = aslinearoperator(AT)
    n = A_linear_operator.shape[0]

    # Initialize the X block with columns of unit 1-norm.
    X = np.ones((n, t))
    if t > 1:
        X[:, 1:] = np.random.randint(0, 2, size=(n, t-1))*2 - 1
    X /= float(n)

    # Iteratively improve the lower bounds.
    # Track extra things, to assert invariants for debugging.
    g_prev = None
    h_prev = None
    k = 1
    ind = range(t)
    while True:
        Y = np.asarray(A_linear_operator.matmat(X))
        g = [norm_1d_1(Y[:, j]) for j in range(t)]
        best_j = np.argmax(g)
        g = sorted(g, reverse=True)
        S = sign_round_up(Y)
        Z = np.asarray(AT_linear_operator.matmat(S))
        h = [norm_1d_inf(row) for row in Z]

        # If this algorithm runs for fewer than two iterations,
        # then its return values do not have the properties indicated
        # in the description of the algorithm.
        # In particular, the entries of g are not 1-norms of any
        # column of A until the second iteration.
        # Therefore we will require the algorithm to run for at least
        # two iterations, even though this requirement is not stated
        # in the description of the algorithm.
        if k >= 2:
            if less_than_or_close(max(h), np.dot(Z[:, best_j], X[:, best_j])):
                break
        h_i_pairs = zip(h, range(n))
        h, ind = zip(*sorted(h_i_pairs, reverse=True)[:t])
        for j in range(t):
            X[:, j] = elementary_vector(n, ind[j])

        # Check invariant (2.2).
        if k >= 2:
            if not less_than_or_close(g_prev[0], h_prev[0]):
                raise Exception('invariant (2.2) is violated')
            if not less_than_or_close(h_prev[0], g[0]):
                raise Exception('invariant (2.2) is violated')

        # Check invariant (2.3).
        if k >= 3:
            for j in range(t):
                if not less_than_or_close(g[j], g_prev[j]):
                    raise Exception('invariant (2.3) is violated')

        # Update for the next iteration.
        g_prev = g
        h_prev = h
        k += 1

    # Return the lower bounds and the corresponding column indices.
    return g, ind


def _onenormest_core(A, AT, t, itmax):
    """
    Compute a lower bound of the 1-norm of a sparse matrix.

    Parameters
    ----------
    A : ndarray or other linear operator
        A linear operator that can produce matrix products.
    AT : ndarray or other linear operator
        The transpose of A.
    t : int, optional
        A positive parameter controlling the tradeoff between
        accuracy versus time and memory usage.
    itmax : int, optional
        Use at most this many iterations.

    Returns
    -------
    est : float
        An underestimate of the 1-norm of the sparse matrix.
    v : ndarray, optional
        The vector such that ||Av||_1 == est*||v||_1.
        It can be thought of as an input to the linear operator
        that gives an output with particularly large norm.
    w : ndarray, optional
        The vector Av which has relatively large 1-norm.
        It can be thought of as an output of the linear operator
        that is relatively large in norm compared to the input.
    nmults : int, optional
        The number of matrix products that were computed.
    nresamples : int, optional
        The number of times a parallel column was observed,
        necessitating a re-randomization of the column.

    Notes
    -----
    This is algorithm 2.4.

    """
    # This function is a more or less direct translation
    # of Algorithm 2.4 from the Higham and Tisseur (2000) paper.
    A_linear_operator = aslinearoperator(A)
    AT_linear_operator = aslinearoperator(AT)
    if itmax < 2:
        raise ValueError('at least two iterations are required')
    if t < 1:
        raise ValueError('at least one column is required')
    n = A.shape[0]
    if t >= n:
        raise ValueError('t should be smaller than the order of A')
    # Track the number of big*small matrix multiplications
    # and the number of resamplings.
    nmults = 0
    nresamples = 0
    # "We now explain our choice of starting matrix.  We take the first
    # column of X to be the vector of 1s [...] This has the advantage that
    # for a matrix with nonnegative elements the algorithm converges
    # with an exact estimate on the second iteration, and such matrices
    # arise in applications [...]"
    X = np.ones((n, t), dtype=float)
    # "The remaining columns are chosen as rand{-1,1},
    # with a check for and correction of parallel columns,
    # exactly as for S in the body of the algorithm."
    if t > 1:
        for i in range(1, t):
            # These are technically initial samples, not resamples,
            # so the resampling count is not incremented.
            resample_column(i, X)
        for i in range(t):
            while column_needs_resampling(i, X):
                resample_column(i, X)
                nresamples += 1
    # "Choose starting matrix X with columns of unit 1-norm."
    X /= float(n)
    # "indices of used unit vectors e_j"
    ind_hist = set()
    est_old = 0
    S = np.zeros((n, t), dtype=float)
    k = 1
    while True:
        Y = np.asarray(A_linear_operator.matmat(X))
        nmults += 1
        mags = [norm_1d_1(Y[:, j]) for j in range(t)]
        est = np.max(mags)
        best_j = np.argmax(mags)
        if est > est_old or k == 2:
            if k >= 2:
                ind_best = ind[best_j]
            w = Y[:, best_j]
        # (1)
        if k >= 2 and est <= est_old:
            est = est_old
            break
        est_old = est
        S_old = S
        if k > itmax:
            break
        S = sign_round_up(Y)
        # (2)
        if every_col_of_X_is_parallel_to_a_col_of_Y(S, S_old):
            break
        if t > 1:
            # "Ensure that no column of S is parallel to another column of S
            # or to a column of S_old by replacing columns of S by rand{-1,1}."
            for i in range(t):
                while column_needs_resampling(i, S, S_old):
                    resample_column(i, S)
                    nresamples += 1
        # (3)
        Z = np.asarray(AT_linear_operator.matmat(S))
        nmults += 1
        h = [norm_1d_inf(row) for row in Z]
        # (4)
        if k >= 2 and max(h) == h[ind_best]:
            break
        # "Sort h so that h_first >= ... >= h_last
        # and re-order ind correspondingly."
        h_i_pairs = zip(h, range(n))
        h, ind = zip(*sorted(h_i_pairs, reverse=True))
        if t > 1:
            # (5)
            # Break if the most promising t vectors have been visited already.
            if set(ind[:t]) <= ind_hist:
                break
            # Put the most promising unvisited vectors at the front of the list
            # and put the visited vectors at the end of the list.
            # Preserve the order of the indices induced by the ordering of h.
            unused_entries = [i for i in ind if i not in ind_hist]
            used_entries = [i for i in ind if i in ind_hist]
            ind = unused_entries + used_entries
        for j in range(t):
            X[:, j] = elementary_vector(n, ind[j])
        ind_hist.update(ind[:t])
        k += 1
    v = elementary_vector(n, ind_best)
    return est, v, w, nmults, nresamples
