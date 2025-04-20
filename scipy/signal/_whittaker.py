import numpy as np
from scipy.linalg import solveh_banded


# TODO:
# 1) C code for _solve_WH_order2_fast
# 2) GCV for lamb
# 3) 2-d, maybe even 3-d WH smoothing

def whittaker_henderson(signal, lamb = 1.0, order=2, weights=None):
    r"""
    Whittaker-Henderson (WH) smoothing/graduation of a discrete signal.

    This implements WH smoothing with a difference penalty of the specified `order` and
    penalty strength `lamb`, see [1] and [2]. WH can be seen as a P-Spline
    (penalized B-Spline) of degree zero for equidistant knots (at the signal
    positions).

    In econometrics, the WH graduation of order 2 is referred to as the Hodrick and
    Prescott filter (https://doi.org/10.2307/2953682).
    
    Parameters
    ----------
    signal : ndarray
        A rank-1 array at least of length 3 representing equidistant data points of a
        signal, e.g. a time series with constant time lag.
    
    lamb : float, optional
        Smoothing or penalty parameter, default is 0.0.

    order : int, default: 2
        The order of the difference penalty, must be at least 1.

    weights : ndarray, option
        A rank-1 array of case weights with the same lenght as `signal`.
        `None` is equivalent to an array of all ones, i.e. `np.ones_like(signal)`.

    Returns
    -------
    x : ndarray
        The WH smoothed signal.

    Notes
    -----
    For the signal :math:`y = (y_1, y_2, \ldots, y_n)` and weights
    :math:`w = (w_1, w_2, \ldots, w_n)`, WH of order :math:`p=2` with smoothing or
    penalty parameter :math:`\lambda` solves the following optimization problem:

    .. math::

        \operatorname{argmin}_{x_i} \sum_i^n w_i (y_i - x_i)^2
        + \lambda \sum_i^{n-p} (\Delta^p x_i) \,,

    with forward difference :math:`\Delta x_i = x_{i+1} - x_i` and
    :math:`\Delta^2 x_i = \Delta(\Delta x_i) = x_{i+2} - 2x_{i+1} + x_i`.
    For every input value :math:`y_i`, it generates a smoothed value :math:`x_i`.

    References
    ----------
    .. [1] Eilers, P.H.C. (2003).
           "A perfect smoother". Analytical Chem. 75, 3631-3636.
           :doi:`10.1021/AC034173T`
    .. [2] Weinert, Howard L. (2007).
           "Efficient computation for Whittaker-Henderson smoothing".
           Computational Statistics and Data Analysis 52:959-74.
           :doi:`10.1016/j.csda.2006.11.038`
    """
    if order < 1 or int(order) != order:
        raise ValueError("Parameter order must be an integer larger equal 1.")

    signal = np.asarray(signal)
    if signal.ndim != 1:
        msg = f"Input signal array must be of shape (n,); got {signal.shape}"
        raise ValueError(msg)

    n = signal.shape[0]
    if n < order + 1:
        msg = f"Input signal array must be at least of shape ({order + 1},); got {n}."
        raise ValueError(msg)

    if weights is not None:
        weights = np.asarray(weights)
        if weights.shape != signal.shape:
            msg = "Parameter weights must have the same shape as the signal array."
            raise ValueError(msg)


    if lamb < 0:
        msg = f"Parameter lamb must be non-negative; got {lamb=}."
        raise ValueError(msg)
    elif lamb == 0.0:
        x = np.asarray(signal).copy()
    else:
        x = _solve_WH_banded(signal, lamb=lamb, order=order, weights=weights)
        # If performance matters and p == 2, think about a C++/Pybind implementation of
        # x = _solve_WH_order2_fast(signal, lamb=lamb)
    return x


def _solve_WH_banded(y, lamb, order=2, weights=None):
    """Solve the WH optimization problem directly with matrices."""
    n = y.shape[0]  # n >= p + 1 was already checked
    p = order  # order of difference penalty
    # Construct penalty matrix M of shape (n-p, n) as if n = 2p+1 (to save memory).
    if n < 2*p + 1:
        M_raw = np.diff(np.eye(n), n=p, axis=0)  # shape (n-p, n)
    else:
        M_raw = np.diff(np.eye(2*p + 1), n=p, axis=0)  # shape (p+1, 2p+1)
    MTM_raw = M_raw.T @ M_raw  # shape (2p+1, 2p+1) if n>=2p+1 else (n, n)
    # Because our matrix A = np.eye(n, dtype=np.float64) + lamb * (M.T @ M) is
    # symmetric and banded with u = l = p, we construct it in the lower "ab"-format
    # for use in solveh_banded, i.e. each row in ab is a subdiagonal of A:
    #   ab[0, :]   = np.diagonal(A, 0)
    #   ab[1, :-1] = np.diagonal(A, 1)
    #   ab[2, :-2] = np.diagonal(A, 2)
    #   ..
    ab = np.zeros((p + 1, min(2*p + 1, n)))
    for i in range(p + 1):
        ab[i, :ab.shape[1] - i] = np.diagonal(MTM_raw, i)
    ab *= lamb
    if n > 2*p + 1:
        ab = np.hstack([
            ab[:, :p+1],
            np.repeat(ab[:, p:p+1], n - (2*p+1), axis=1),
            ab[:, -p:],
        ])
    if weights is None:
        ab[0, :] += 1.0  # This corresponds to np.eye(n).
        x = solveh_banded(ab, y, lower=True)
    else:
        ab[0, :] += weights
        x = solveh_banded(ab, weights * y, lower=True)
    return x


def _solve_WH_order2_fast(y, lamb):
    """Efficiently solve WH of order 2 according to Weinert.

    Needs order = 2 and n >= 3 (data points).

    Weinert (2007)
    "Efficient computation for Whittaker-Henderson smoothing".
    Computational Statistics and Data Analysis 52:959-74.
    https://doi.org/10.1016/j.csda.2006.11.038
    """
    n = y.shape[0]
    # Convert penalty to convention of Weinert (2007), i.e. A = lambda I + M'M.
    lamb = 1 / lamb
    # A = LDL' decompositon, i.e. L is unit lower triangular and banded
    # (bandwith=2) and D diagonal.
    # First subdiagonal of L is (-e_1, .., -e_{n-1}) and 2nd subdiagonal
    # (f_1, .., f_{n-2}). Diagonal of D is (d_1, .., d_n).
    # The equation A @ x = lamb * y becomes
    # LD @ b = lamb * y  (I)
    # L.T @ x = b        (II)
    # We shift Weinert's 1-based indices to 0-based indices, Eq. 2.2-2.6
    # Solve problem (I)
    b = np.empty(n)
    e = np.empty(n)
    f = np.empty(n)
    # i=0
    d = 1 + lamb
    f[0] = 1 / d
    mu = 2
    e[0] = mu * f[0]
    b[0] = f[0] * lamb * y[0]
    mu_old = mu
    # i=1
    if n == 3:
        d = 4 + lamb - mu_old * e[0]
        mu = 2 - e[0]
    else:
        d = 5 + lamb - mu_old * e[0]
        mu = 4 - e[0]
    f[1] = 1 / d
    e[1] = mu * f[1]
    b[1] = f[1] * (lamb * y[1] + mu_old * b[0])
    mu_old = mu
    for i in range(2, n-2):
        d = 6 + lamb - mu_old * e[i-1] - f[i-2]
        f[i] = 1 / d
        mu = 4 - e[i-1]
        e[i] = mu * f[i]
        b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])
        mu_old = mu
    # i=n-2
    if n >= 4:
        i = n - 2
        d = 5 + lamb - mu_old * e[i-1] - f[i-2]
        f[i] = 1 / d
        mu = 2 - e[i-1]
        e[i] = mu * f[i]
        b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])
        mu_old = mu
    # i=n-1
    i = n - 1
    d = 1 + lamb - mu_old * e[i-1] - f[i-2]
    f[i] = 1 / d
    b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])

    # Solve problem (II)
    x = np.empty(n)
    x[n-1] = b[n-1]
    x[n-2] = b[n-2] + e[n-2] * x[n-1]
    for i in range(n-3, -1, -1):
        x[i] = b[i] + e[i] * x[i+1] - f[i] * x[i+2]
    return x
