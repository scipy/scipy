import numpy as np
from scipy.linalg import solve, toeplitz


# TODO:
# 1) Add parameter `order` for the order of the difference penalty
# 2) C code for _solve_WH_order2_fast
# 3) Use sparse/banded arrays for _solve_WH_order_direct
# 4) GCV for lamb
# 5) Case weights
# 6) 2-d, maybe even 3-d WH smoothing

def whittaker_henderson(signal, lamb = 0.0):
    r"""
    Whittaker-Henderson (WH) smoothing/graduation of a discrete signal.

    This implements WH of order 2, see [1] and [2]. WH can be seen as a P-Spline
    (penalized B-Spline) of degree zero for equidistant knots.

    In econometrics, the WH graduation of order 2 is referred to as the Hodrick and
    Prescott filter (https://doi.org/10.2307/2953682).
    
    Parameters
    ----------
    signal : ndarray
        A rank-1 array at least of length 3 representing equidistant data points of a
        signal, e.g. a time series with constant time lag.
    
    lamb : float, optional
        Smoothing or penalty parameter, default is 0.0.

    Returns
    -------
    x : ndarray
        WH smoothed signal.

    Notes
    -----
    For signal of :math:`y = (y_1, y_2, \ldots, y_n)`, WH of order :math:`p=2` with
    smoothing parameter :math:`\lambda` solves the following optimization problem:

    .. math::

        \operatorname{argmin}_{x_i} \sum_i^n w_i (y_i - x_i)^2
        + \lambda \sum_i^{n-p} (\Delta^p x_i) \,,

    with :math:`\Delta x_i = x_{i+1} - x_i` and
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
    signal = np.asarray(signal)
    if signal.ndim != 1:
        msg = f"Input signal array must be of shape (n,); got {signal.shape}"
        raise ValueError(msg)

    n = signal.shape[0]
    if n < 3:
        # signal must be at least of length 3 (=order+1).
        msg = f"Input signal array must be at least of shape (3,); got {n}"
        raise ValueError(msg)

    if lamb < 0:
        msg = f"Parameter lamb must be non-negative; got {lamb=}."
        raise ValueError(msg)
    elif lamb == 0.0:
        x = np.asarray(signal).copy()
    else:
        if n < 5:
            x = _solve_WH_order_direct(signal, lamb=lamb)
        else:
            x = _solve_WH_order2_fast(signal, lamb=lamb)
    return x


def _solve_WH_order_direct(y, lamb):
    """Solve the WH optimization problem directly with matrices."""
    n = y.shape[0]
    p = 2  # order of difference penalty
    col = np.zeros(n - p)
    col[0] = 1
    row = np.zeros(n)
    row[:3] = [1, -2, 1]
    M = toeplitz(c=col, r=row)
    A = np.eye(n, dtype=np.float64) + lamb * (M.T @ M)
    x = solve(A, y)
    return x


def _solve_WH_order2_fast(y, lamb):
    """Efficiently solve WH of order 2 according to Weinert.

    Needs order = 2 and n >= 5 (data points).

    Weinert (2007)
    "Efficient computation for Whittaker-Henderson smoothing".
    Computational Statistics and Data Analysis 52:959-74.
    https://doi.org/10.1016/j.csda.2006.11.038
    """
    n = y.shape[0]
    # Convert penalty to convention of Weinert (2007)
    lamb = 1 / lamb
    # A = LDL decompositon, i.e. L is unit lower triangular and banded
    # (bandwith=2) and D diagonal.
    # The equation A @ x = lamb * y becomes
    # LD @ b = lamb * y  (I)
    # L.T @ x = b        (II)
    # We shift Weinert's 1-based indices to 0-based indices, Eq. 2.2-2.6
    # problem (I)
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
    d = 5 + lamb - mu_old * e[0]
    f[1] = 1 / d
    mu = 4 - e[0]
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

    # problem (II)
    x = np.empty(n)
    x[n-1] = b[n-1]
    x[n-2] = b[n-2] + e[n-2] * x[n-1]
    for i in range(n-3, -1, -1):
        x[i] = b[i] + e[i] * x[i+1] - f[i] * x[i+2]
    return x
