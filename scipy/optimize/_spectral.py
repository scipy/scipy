"""
Spectral Algorithm for Nonlinear Equations
"""
from __future__ import division, absolute_import, print_function

import collections

import numpy as np
from scipy.optimize import OptimizeResult
from scipy.optimize.optimize import _check_unknown_options
from .linesearch import _nonmonotone_line_search_cruz, _nonmonotone_line_search_cheng

class _NoConvergence(Exception):
    pass


def _root_df_sane(func, x0, args=(), ftol=1e-8, fatol=1e-300, maxfev=1000,
                  fnorm=None, callback=None, disp=False, M=10, eta_strategy=None,
                  sigma_eps=1e-10, sigma_0=1.0, line_search='cruz',
                  stagnation_limit=100, **unknown_options):
    r"""
    Solve nonlinear equation with the DF-SANE method

    Parameters
    ----------
    func : callable
        Function whose root to find
    x0 : ndarray
        Initial guess
    args : tuple
        Extra parameters to `func`
    ftol : float, optional
        Relative norm tolerance.
    fatol : float, optional
        Absolute norm tolerance.
        Algorithm terminates when ``||func(x)|| < fatol + ftol ||func(x_0)||``.
    fnorm : callable, optional
        Norm to use in the convergence check. If None, 2-norm is used.
    callback : callable, optional
        Callback function to call on each iteration. Called as
        ``callable(x, f)`` where ``x`` is the current iterate and
        ``f`` the residual.
    maxfev : int, optional
        Maximum number of function evaluations.
    disp : bool, optional
        Whether to print convergence process to stdout.
    eta_strategy : callable, optional
        Choice of the ``eta_k`` parameter, which gives slack for growth
        of ``||F||_2**2``.  Called as ``eta_k = eta_strategy(k, x, F)`` with
        `k` the iteration number, `x` the current iterate and `F` the current
        residual. Should satisfy ``eta_k > 0`` and ``sum(eta, k=0..inf) < inf``.
        Default: ``||F_0||_2**2 / (1 + k)**2`` where F_0 is the initial residual.
    sigma_eps : float, optional
        The spectral coefficient is constrained to ``sigma_eps < sigma < 1/sigma_eps``.
        Default: 1e-10
    sigma_0 : float, optional
        Initial spectral coefficient.
        Default: 1.0
    M : int, optional
        Number of iterates to include in the 'cruz' nonmonotonic line search.
        Default: 10
    stagnation_limit : int, optional
        If merit function does not decrease in `stagnation_limit` iterations,
        terminate with failure.
        Default: 100
    line_search : {'cruz', 'cheng'}, optional
        Type of line search to employ. 'cruz' is the original one defined in [1],
        'cheng' is a modified search defined in [3].
        Default: 'cruz'

    References
    ----------
    .. [1] "Spectral residual method without gradient information for solving
           large-scale nonlinear systems of equations." W. La Cruz,
           J.M. Martinez, M. Raydan. Math. Comp. **75**, 1429 (2006).
    .. [2] W. La Cruz, Opt. Meth. Software, 29, 24 (2014).
    .. [3] W. Cheng, D.-H. Li. IMA J. Numer. Anal. **29**, 814 (2009).

    """
    _check_unknown_options(unknown_options)

    if line_search not in ('cheng', 'cruz'):
        raise ValueError("Invalid value %r for 'line_search'" % (line_search,))

    nexp = 2

    if eta_strategy is None:
        # Different choice from [1], as their eta is not invariant
        # vs. scaling of F.
        def eta_strategy(k, x, F):
            # Obtain squared 2-norm of the initial residual from the outer scope
            return f_0 / (1 + k)**2

    if fnorm is None:
        def fnorm(F):
            # Obtain squared 2-norm of the current residual from the outer scope
            return f_k**(1.0/nexp)

    def fmerit(x, F):
        return np.linalg.norm(F.ravel())**nexp

    nfev = [0]
    f, x_k, x_shape, f_k, F_k, is_complex = _wrap_func(func, x0, fmerit, nfev, maxfev, args)

    k = 0
    f_0 = f_k
    sigma_k = sigma_0

    F_0_norm = fnorm(F_k)

    # For the 'cruz' line search
    prev_fs = collections.deque([f_k], M)

    # For the 'cheng' line search
    Q = 1.0
    C = f_0

    # Track improvement
    last_improved_k = 0
    last_improved_f = f_0

    converged = False
    message = "too many function evaluations required"

    while True:
        F_k_norm = fnorm(F_k)

        # Control spectral parameter, from [2]
        if abs(sigma_k) > 1/sigma_eps:
            sigma_k = 1/sigma_eps * np.sign(sigma_k)
        elif abs(sigma_k) < sigma_eps:
            sigma_k = sigma_eps

        if disp:
            print("iter %d: ||F|| = %g, sigma = %g" % (k, F_k_norm, sigma_k))

        if callback is not None:
            callback(x_k, F_k)

        # Check convergence
        if F_k_norm < ftol * F_0_norm + fatol:
            # Converged!
            message = "successful convergence"
            converged = True
            break

        if f_k < last_improved_f:
            last_improved_f = f_k
            last_improved_k = k
        elif k - last_improved_k > stagnation_limit:
            message = "convergence stagnated"
            converged = False
            break

        # Line search direction
        d = -sigma_k * F_k

        # Nonmonotone line search
        eta = eta_strategy(k, x_k, F_k)
        try:
            if line_search == 'cruz':
                alpha, xp, fp, Fp = _nonmonotone_line_search_cruz(f, x_k, d, prev_fs, eta=eta)
            elif line_search == 'cheng':
                alpha, xp, fp, Fp, C, Q = _nonmonotone_line_search_cheng(f, x_k, d, f_k, C, Q, eta=eta)
        except _NoConvergence:
            break

        # Update spectral parameter
        s_k = xp - x_k
        y_k = Fp - F_k
        sigma_k = np.vdot(s_k, s_k) / np.vdot(s_k, y_k)

        # Take step
        x_k = xp
        F_k = Fp
        f_k = fp

        # Store function value
        if line_search == 'cruz':
            prev_fs.append(fp)

        k += 1

    x = _wrap_result(x_k, is_complex, shape=x_shape)
    F = _wrap_result(F_k, is_complex)

    result = OptimizeResult(x=x, success=converged,
                            message=message,
                            fun=F, nfev=nfev[0], nit=k)

    return result


def _fixpoint_squarem(func, x0, args=(), xtol=1e-8, xatol=1e-300,
                      maxfev=1000, fnorm=None, callback=None, disp=False,
                      fmerit=None, step_rule=3, **unknown_options):
    r"""
    Solve a fixed-point equation with the SQUAREM method.

    Parameters
    ----------
    func : callable
        Function whose fixed point to find
    x0 : ndarray
        Initial guess
    args : tuple
        Extra parameters to `func`
    callback : callable, optional
        Callback function to call on each iteration. Called as
        ``callable(x, f)`` where ``x`` is the current iterate and
        ``f`` the residual.
    xtol : float, optional
        Relative norm tolerance.
    xatol : float, optional
        Absolute norm tolerance.
        Algorithm terminates when ``||func(x) - x|| < xtol ||x|| + xatol``.
    fnorm : callable, optional
        Norm to use in the convergence check. If None, 2-norm is used.
    maxfev : int, optional
        Maximum number of function evaluations.
    disp : bool, optional
        Whether to print convergence process to stdout.
    step_rule : {1, 2, 3}, optional
        Which of the SQUAREM step length rules to use, see below.
        Default: 2
    fmerit : callable, optional
        Merit function to be minimized in the globalization step;
        called as ``fmerit(x, F)``. For EM, take the negative of
        the likelihood function. Default: ||F - x||**2

    Notes
    -----
    The convergence result for the algorithm in [1] assumes that the function
    is a contraction in the norm specified by the merit function.

    There are three possible choices for the step length, see [1]_::

        (1)  alpha_k = vdot(r, v) / vdot(v, v)
        (2)  alpha_k = vdot(r, r) / vdot(r, v)
        (3)  alpha_k = -sqrt(vdot(r, r) / vdot(v, v))

    The different choices minimize different step length model objective functions.

    References
    ----------
    .. [1] R. Varadhan, C. Roland. Scand. J. Statistics, 35, 335 (2008).

    """
    _check_unknown_options(unknown_options)

    if step_rule not in (1, 2, 3):
        raise ValueError("Invalid step_rule %r; shoul be one of (1, 2, 3)" % (step_rule,))

    if fmerit is None:
        def fmerit(x, F):
            return np.linalg.norm((x - F).ravel())**2

    nfev = [0]
    f, x_k, x_shape, f_k, F_k, is_complex = _wrap_func(func, x0, fmerit, nfev, maxfev, args)

    k = 0
    alpha_k = 1.0
    backtrack_theta = 0.1

    if fnorm is None:
        def fnorm(z):
            return np.linalg.norm(z)

    converged = False
    message = "Too many function evaluations required"

    while True:
        dx_norm = fnorm(F_k - x_k)

        try:
            _, FF_k = f(F_k)
        except _NoConvergence:
            break

        # SQUAREM spectral parameter
        r_k = F_k - x_k
        v_k = (FF_k - F_k) - r_k

        with np.errstate(invalid='ignore', over='ignore'):
            if step_rule == 1:
                alpha_k = np.vdot(r_k, v_k) / np.vdot(v_k, v_k)
            elif step_rule == 2:
                alpha_k = np.vdot(r_k, r_k) / np.vdot(r_k, v_k)
            elif step_rule == 3:
                alpha_k = -np.linalg.norm(r_k) / np.linalg.norm(v_k)

        if disp:
            print("iter %d: ||dx|| = %g, sigma = %g" % (k, dx_norm, alpha_k))

        if callback is not None:
            callback(x_k)

        if dx_norm < xtol * fnorm(x_k) + xatol:
            # Converged!
            message = "Converged successfully"
            converged = True
            break

        # Backtracking
        try:
            # The approach from [1]
            if not (alpha_k <= -1):
                alpha_k = -1

            while True:
                xp = x_k - 2 * alpha_k * r_k + alpha_k**2 * v_k
                fp, Fp = f(xp)
                if fp < f_k or not (alpha_k < -1):
                    break
                alpha_k = backtrack_theta * alpha_k + (-1)*(1 - backtrack_theta)
        except _NoConvergence:
            break

        # Stabilization
        xp = Fp

        # Take step
        x_k = xp
        f_k, F_k = f(xp)

        k += 1

    x = _wrap_result(x_k, is_complex, shape=x_shape)
    F = _wrap_result(F_k, is_complex)

    result = OptimizeResult(x=x, success=converged,
                            message=message,
                            fun=F, nfev=nfev[0], nit=k)

    return result


def _wrap_func(func, x0, fmerit, nfev_list, maxfev, args=()):
    """
    Wrap a function and an initial value so that (i) complex values
    are wrapped to reals, and (ii) value for a merit function
    fmerit(x, f) is computed at the same time, (iii) iteration count
    is maintained and an exception is raised if it is exceeded.

    Parameters
    ----------
    func : callable
        Function to wrap
    x0 : ndarray
        Initial value
    fmerit : callable
        Merit function fmerit(f) for computing merit value from residual.
    nfev_list : list
        List to store number of evaluations in. Should be [0] in the beginning.
    maxfev : int
        Maximum number of evaluations before _NoConvergence is raised.
    args : tuple
        Extra arguments to func

    Returns
    -------
    wrap_func : callable
        Wrapped function, to be called as
        ``F, fp = wrap_func(x0)``
    x0_wrap : ndarray of float
        Wrapped initial value; raveled to 1D and complex
        values mapped to reals.
    x0_shape : tuple
        Shape of the initial value array
    f : float
        Merit function at F
    F : ndarray of float
        Residual at x0_wrap
    is_complex : bool
        Whether complex values were mapped to reals

    """
    x0 = np.asarray(x0)
    x0_shape = x0.shape
    F = np.asarray(func(x0, *args)).ravel()
    is_complex = np.iscomplexobj(x0) or np.iscomplexobj(F)
    x0 = x0.ravel()

    nfev_list[0] = 1

    if is_complex:
        def wrap_func(x):
            if nfev_list[0] >= maxfev:
                raise _NoConvergence()
            nfev_list[0] += 1
            z = _real2complex(x).reshape(x0_shape)
            v = np.asarray(func(z, *args)).ravel()
            f = fmerit(z, v)
            F = _complex2real(v)
            return f, F

        x0 = _complex2real(x0)
        F = _complex2real(F)
    else:
        def wrap_func(x):
            if nfev_list[0] >= maxfev:
                raise _NoConvergence()
            nfev_list[0] += 1
            z = x.reshape(x0_shape)
            F = np.asarray(func(z, *args))
            f = fmerit(z, F)
            F = F.ravel()
            return f, F

    return wrap_func, x0, x0_shape, fmerit(x0, F), F, is_complex


def _wrap_result(result, is_complex, shape=None):
    """
    Convert from real to complex and reshape result arrays.
    """
    if is_complex:
        z = _real2complex(result)
    else:
        z = result
    if shape is not None:
        z = z.reshape(shape)
    return z


def _real2complex(x):
    return np.ascontiguousarray(x, dtype=float).view(np.complex128)


def _complex2real(z):
    return np.ascontiguousarray(z, dtype=complex).view(np.float64)
