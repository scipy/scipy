import numpy as np
from scipy.linalg import LinAlgError
from scipy.linalg.lapack import get_lapack_funcs
from scipy.optimize import minimize_scalar
from scipy.special import binom

# TODO:
# 1) C code for _solve_WH_order2_fast
# 2) GCV for lamb
# 3) 2-d, maybe even 3-d WH smoothing


def _solveh_banded(ab, b, calc_logdet=False):
    """
    Solve the equation ``a @ x = b`` for ``x``,  where ``a`` is the 
    Hermitian positive-definite banded matrix defined by `ab`.

    Same as scipy.linalg.solveh_banded(lower=True, check_finite=False), but:
    - also returns the log of the determinant
    - no input validation
    - only real values, no complex
    - only `lower = True` code path
    - always overwrite_XX = False
    - b only a 1-dim array

    Parameters
    ----------
    ab : (``u`` + 1, M) array_like
        Banded matrix
    b : (M,) array_like
        Right-hand side

    Returns
    -------
    x : (M,) ndarray
        The solution to the system ``a x = b``. Shape of return matches shape of `b`.
    logdet : float
        Logarithm of the determinant of `ab`.
    """
    a1 = ab
    b1 = b
    overwrite_b = False
    overwrite_ab = False
    logdet = 0

    if a1.shape[0] == 2:
        ptsv, = get_lapack_funcs(("ptsv",), (a1, b1))
        # We assume lower=True and real arrays
        d = a1[0, :]
        e = a1[1, :-1]
        # ptsv uses LDL', returnes d=diag(D), du=diag(L, -1)
        d, du, x, info = ptsv(d, e, b1, overwrite_ab, overwrite_ab, overwrite_b)
        if calc_logdet:
            logdet = np.log(d).sum()
    else:
        pbsv, = get_lapack_funcs(("pbsv",), (a1, b1))
        # pbsv uses Cholesky LL', returns c=L in ab-storage format
        c, x, info = pbsv(a1, b1, lower=True, overwrite_ab=overwrite_ab,
                          overwrite_b=overwrite_b)
        if calc_logdet:
            logdet = 2 * np.log(c[0, :]).sum()
    if info > 0:
        raise LinAlgError(f"{info}th leading minor not positive definite")
    if info < 0:
        raise ValueError(f"illegal value in {-info}th argument of internal pbsv")
    return x, logdet


def whittaker_henderson(signal, lamb="reml", order=2, weights=None):
    r"""
    Whittaker-Henderson (WH) smoothing/graduation of a discrete signal.

    This implements WH smoothing with a difference penalty of the specified `order` and
    penalty strength `lamb`, see [1]_, [2]_ and [3]_. WH can be seen as a P-Spline
    (penalized B-Spline) of degree zero for equidistant knots (at the signal
    positions).

    In econometrics, the WH graduation of order 2 is referred to as the Hodrick and
    Prescott filter [4]_.
    
    Parameters
    ----------
    signal : ndarray
        A rank-1 array at least of length `order + 1` representing equidistant data
        points of a signal, e.g. a time series with constant time lag.
    
    lamb : str or float, optional
        Smoothing or penalty parameter, default is "reml" which minimized the REML
        criterion to find the parameter `lamb`. If a number is passed, it must be
        non-negative and it is used directly.

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
    :math:`w = (w_1, w_2, \ldots, w_n)`, WH of order :math:`p` with smoothing or
    penalty parameter :math:`\lambda` solves the following optimization problem:

    .. math::

        \operatorname{argmin}_{x_i} \sum_i^n w_i (y_i - x_i)^2
        + \lambda \sum_i^{n-p} (\Delta^p x_i)^2 \,,

    with forward difference :math:`\Delta x_i = x_{i+1} - x_i` and
    :math:`\Delta^2 x_i = \Delta(\Delta x_i) = x_{i+2} - 2x_{i+1} + x_i`.
    For every input value :math:`y_i`, it generates a smoothed value :math:`x_i`.

    References
    ----------
    .. [1] Whittaker-Henderson smoothing,
           https://en.wikipedia.org/wiki/Whittaker%E2%80%93Henderson_smoothing
    .. [2] Eilers, P.H.C. (2003).
           "A perfect smoother". Analytical Chem. 75, 3631-3636.
           :doi:`10.1021/AC034173T`
    .. [3] Weinert, Howard L. (2007).
           "Efficient computation for Whittaker-Henderson smoothing".
           Computational Statistics and Data Analysis 52:959-74.
           :doi:`10.1016/j.csda.2006.11.038`
    .. [4] Hodrick, R. J., & Prescott, E. C. (1997).
           Postwar U.S. Business Cycles: An Empirical Investigation.
          :doi:`10.2307/2953682`
    """
    if order < 1 or int(order) != order:
        raise ValueError("Parameter order must be an integer larger equal 1.")

    signal = np.asarray(signal)
    if signal.ndim != 1:
        msg = f"Input array signal must be of shape (n,); got {signal.shape}"
        raise ValueError(msg)

    n = signal.shape[0]
    if n < order + 1:
        msg = f"Input array signal must be at least of shape ({order + 1},); got {n}."
        raise ValueError(msg)

    if weights is not None:
        weights = np.asarray(weights)
        if weights.shape != signal.shape:
            msg = "Input array weights must have the same shape as the signal array."
            raise ValueError(msg)
        
    if weights is None:
        if not np.isfinite(signal).all():
            raise ValueError("Input array signal must be finite.")
    else:
        if not np.isfinite(weights).all():
            raise ValueError("Input array weights must be finite.")
        if (mask := ~np.isfinite(signal)).any():
            # Only weights * y matter in the end and weights must be 0 for all
            # non-finite elements of signal.
            if not (weights[mask] == 0).all():
                raise ValueError("Input array weights must be zero for all non-finite "
                                 "elements of signal.")
            signal = np.nan_to_num(signal)


    msg = f"Parameter lamb must be string 'reml' or a non-negative float; got {lamb=}."
    if isinstance(lamb, str):
        if lamb != "reml":
            raise ValueError(msg)
        def criterion(loglamb):
            return -_reml(lamb=np.exp(loglamb), y=signal, order=order, weights=weights)
        opt = minimize_scalar(criterion, bracket=[-10, 10])
        lamb = np.exp(opt.x)

    if lamb < 0:
        raise ValueError(msg)
    elif lamb == 0.0:
        x = np.asarray(signal).copy()
    else:
        # If performance matters and p == 2, think about a C++/Pybind implementation of
        # x = _solve_WH_order2_fast(signal, lamb=lamb)
        x, _ = _solve_WH_banded(signal, lamb=lamb, order=order, weights=weights)
    return x


def _solve_WH_banded(y, lamb, order=2, weights=None, calc_logdet=False):
    """
    Solve the WH optimization problem via the normal equations.
    
    A @ x = y
    A = I + lamb * P = I + lamb * D' @ D
    D = difference matrix of order=`order` 

    With weights W = diag(weights):
    A = W + lamb * P
    A @ x = W @ y

    Returns
    -------
    x : ndarray
        The solution.
    logdet : float
        Logarithm of the determinant of matrix A.
    """
    n = y.shape[0]  # n >= p + 1 was already checked
    p = order  # order of difference penalty
    # Construct penalty matrix P = D'D of shape (n-p, n) as if n = 2p+1 (to save
    # memory).
    if n < 2*p + 1:
        D = np.diff(np.eye(n), n=p, axis=0)  # shape (n-p, n)
    else:
        D = np.diff(np.eye(2*p + 1), n=p, axis=0)  # shape (p+1, 2p+1)
    P_raw = D.T @ D  # shape (2p+1, 2p+1) if n >= 2p+1 else (n, n)
    # Because our matrix A = np.eye(n, dtype=np.float64) + lamb * (D.T @ D) is
    # symmetric and banded with u = l = p, we construct it in the lower "ab"-format
    # for use in solveh_banded, i.e. each row in ab is a subdiagonal of A:
    #   ab[0, :]   = np.diagonal(A, 0)
    #   ab[1, :-1] = np.diagonal(A, 1)
    #   ab[2, :-2] = np.diagonal(A, 2)
    #   ..
    ab = np.zeros((p + 1, min(2*p + 1, n)))
    for i in range(p + 1):
        ab[i, :ab.shape[1] - i] = np.diagonal(P_raw, i)
    ab *= lamb
    if n > 2*p + 1:
        ab = np.hstack([
            ab[:, :p+1],
            np.repeat(ab[:, p:p+1], n - (2*p+1), axis=1),
            ab[:, -p:],
        ])
    if weights is None:
        ab[0, :] += 1.0  # This corresponds to np.eye(n).
        x, logdet = _solveh_banded(ab, y, calc_logdet=calc_logdet)
    else:
        ab[0, :] += weights
        x, logdet = _solveh_banded(ab, weights * y, calc_logdet=calc_logdet)
    return x, logdet


def _solve_WH_order2_fast(y, lamb):
    """Efficiently solve WH of order 2 according to Weinert.

    Needs order = 2 and n >= 3 (data points).

    Weinert (2007)
    "Efficient computation for Whittaker-Henderson smoothing".
    Computational Statistics and Data Analysis 52:959-74.
    https://doi.org/10.1016/j.csda.2006.11.038
    """
    n = y.shape[0]
    # Convert penalty to convention of Weinert (2007), i.e. A = lambda I + D'D.
    # Note that Weinert denotes the difference matrix M instead of D.
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


def _logdet_difference_matrix(order, n):
    """Logarithm of the determinant of the difference matrix.

    If D is the difference matrix of order=p, then this computes `log det(D @ D.T)`
    which equals the log of the sum of non-zero eigenvalues of `D.T @ D`.
    """
    # product of eigenvalues =
    # prod(binom(n+i-1, 2i-1), i=1..p) / prod(binom(2i, i), i=1..p-1)
    # How to derive this formula? Well, ... some magic.
    p = order
    if order == 1:
        return np.log(n)
    logdet = 0.0
    for i in range(1, p+1):
        logdet += np.log(binom(n + i - 1, 2*i - 1) / binom(2*i, i)) 
    logdet += np.log(binom(2*p, p))
    return logdet


def _reml(lamb, y, order, weights=None):
    """Calculate the restictricted maximum likelihood (REML).
    
    Parameters
    ----------
    lamb : penalty

    y : signal

    x : smoothed signal

    order : oder of the difference penalty.

    weights : case weights

    Returns
    -------
    reml : REML criterion

    References
    ----------
    - Biessy https://arxiv.org/abs/2306.06932 (version 4)
    - Wood https://doi.org/10.1111/j.1467-9868.2010.00749.x
    """
    n = y.shape[0]
    x, logdet = _solve_WH_banded(
        y=y, lamb=lamb, order=order, weights=weights, calc_logdet=True
    )
    logdet_DtD = _logdet_difference_matrix(order=order, n=n)
    residual = y - x
    # Eq. 12 of Biessy gives the REML criterion:
    # REML(lambda, sigma) = (log of restriced maximum likelihood)
    #     = -1/2 ((y - theta) W (y - theta) / sigma^2 + lambda theta D'D theta / sigma^2
    #             - log|lambda D'D| + log|(W + lambda D'D)| + (n - p) log(sigma^2)
    #             + const
    #            )
    # where the constant term "const" does not depend on lambda or sigma and p is the
    # order of the difference penalty.
    # Note that Biessy then does not mention to use the profiled REML criterion, i.e.,
    # analytically plug in the optimal sigma^2. This gives us
    #     sigma^2 = r2 / (n - p)
    #     r^2     = (y - theta) W (y - theta) + lambda theta D'D theta
    #     profiled REML(lambda) =
    #         -1/2 (
    #               (n-p) (1 + log(r^2 / (n-p)))
    #               -log|lambda D'D| + log|W + lambda D'D| + const
    #              )
    # This can be compared to Eq. 41 of Bates et al
    # https://doi.org/10.18637/jss.v067.i01.
    # An alternative derivation stems from a mixed model formulation of P-splines, see
    # Currie and Durban https://doi.org/10.1191/1471082x02st039ob or Boer
    # https://doi.org/10.1177/1471082X231178591. One then has 2 variance parameters
    # sigma^2 (from y) and tau^2 (from the random effect) leading to
    #     -1/2 ((y - theta) W (y - theta) / sigma^2 + theta D'D theta / tau^2 + ...
    # One then sets tau^2 = sigma^2 / lambda.
    if weights is None:
        r2 = residual @ residual
    else:
        r2 = residual @ (weights * residual)
    r2 += lamb * np.sum(np.diff(x, n=order)**2)  # + lambda theta D'D theta
    reml = (n - order) * (1 + np.log(r2 / (n - order)))
    reml -= (n - order) * np.log(lamb) + logdet_DtD  # -log|lambda D'D|
    reml += logdet  # +log|W + lambda D'D|
    reml *= -0.5
    return reml
