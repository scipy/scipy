"""Sparse block 1-norm estimator.
"""

from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
import scipy.sparse.linalg

__all__ = ['condest']



def condest(A, t=2, itmax=5, compute_v=False, compute_w=False):
    """
    Compute a lower bound of the 1-norm of a sparse matrix.

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
    nmults : int, optional
        The number of matrix products that were computed.
    nresamples : int, optional
        The number of times a parallel column was observed,
        necessitating a re-randomization of the column.

    Notes
    -----
    This is algorithm 2.4 of [1]_.

    In [2]_ it is described as follows.
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
    est, v, w, nmults, nresamples = _condest_core(A, A.T, t, itmax)
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
    #return np.array_equal(v, w) or np.array_equal(v, -w)
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

def _condest_core(A, AT, t, itmax):
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
        Larger values take longer and use more memory
        but give more accurate output.
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

    """
    # This function is a more or less direct translation
    # of Algorithm 2.4 from the Higham and Tisseur (2000) paper.
    A_linear_operator = scipy.sparse.linalg.aslinearoperator(A)
    AT_linear_operator = scipy.sparse.linalg.aslinearoperator(AT)
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
    ind = [0]*n
    S = np.zeros((n, t), dtype=float)
    k = 1
    while True:
        Y = A_linear_operator.matmat(X)
        nmults += 1
        mags = [norm_1d_1(Y[:, j]) for j in range(t)]
        est = np.max(mags)
        best_j = np.argmax(mags)
        if est > est_old or k == 2:
            ind_best = ind[best_j]
            # is this wrong in the paper?
            #w = Y[:, ind_best]
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
        Z = AT_linear_operator.matmat(S)
        nmults += 1
        h = [norm_1d_inf(row) for row in Z]
        # (4)
        if k >= 2 and max(h) == h[ind_best]:
            break
        # "Sort h so that h_first >= ... >= h_last
        # and re-order ind correspondingly."
        h, ind = zip(*sorted(zip(h, range(n)), reverse=True))
        if t > 1:
            # (5)
            # break if ind is contained in ind_hist
            if set(ind) <= ind_hist:
                break
            # Replace the first t entries of ind
            # by the first t indices in ind that are not in ind_hist.
            # ...
            # What if all but fewer than t entries of ind are in ind_hist?
            # Is ind supposed to be a permutation?
            # If so, then why is it initialized to all zeros?
            unused_entries = []
            for entry in ind:
                if entry not in ind_hist:
                    unused_entries.append(entry)
                    if len(unused_entries) == t:
                        break
            if len(unused_entries) != t:
                # This is a corner case that was not addressed in the paper
                # and which I will treat like condition (5).
                break
            # If ind is supposed to be a permutation,
            # which I cannot discern from the paper,
            # then this should "swap" instead of "replace."
            ind = unused_entries + list(ind[t:])
        for j in range(t):
            X[:, j] = elementary_vector(n, ind[j])
        ind_hist.update(ind[:t])
        k += 1
    v = elementary_vector(n, ind_best)
    return est, v, w, nmults, nresamples


