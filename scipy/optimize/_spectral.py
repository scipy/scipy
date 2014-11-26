"""
Spectral Algorithm for Nonlinear Equations
"""
from __future__ import division, absolute_import, print_function

import warnings
import numpy as np
from scipy.optimize import OptimizeResult
from scipy.optimize.optimize import _check_unknown_options, OptimizeWarning
from .linesearch import _nonmonotone_line_search_cruz

class _NoConvergence(Exception):
    pass


def _root_df_sane(func, x0, args=(), jac=None, ftol=1e-8, fatol=1e-300, maxfev=1000,
                  fnorm=None, callback=None, disp=False, M=None, eta_strategy=None,
                  **unknown_options):
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
    jac
        Unused, this is a derivative-free solver
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
    eta_stragegy : callable, optional
        Choice of the ``eta_k`` parameter, which gives slack for growth
        of ``||F||**2``.  Called as ``eta_k = eta_strategy(k, x, F)`` with
        `k` the iteration number, `x` the current iterate and `F` the current
        residual. Should satisfy ``eta_k > 0`` and ``sum(eta, k=0..inf) < inf``.
        Default: ``||F||**2 / (1 + k)**2``.
    M : int, optional
        Number of iterates to include in the nonmonotonic line search.
        Default: 10

    References
    ----------
    .. [1] "Spectral residual method without gradient information for solving
           large-scale nonlinear systems of equations." W. La Cruz,
           J.M. Martinez, M. Raydan. Math. Comp. **75**, 1429 (2006).
    .. [2] W. La Cruz, Opt. Meth. Software, 29, 24 (2014).

    """
    _check_unknown_options(unknown_options)
    if jac is not None:
        warnings.warn("DF-SANE solver does not use a Jacobian", category=OptimizeWarning)

    nexp = 2
    sigma_eps = 1e-10

    if M is None:
        M = 10

    if eta_strategy is None:
        # Different choice from [1], as their eta is not invariant
        # vs. scaling of F.
        #
        # This choice satisfies the requirements in the paper.  f_k is
        # bounded: the maximum allowed by the line search is
        #
        #     f_{k+1} = f_k * [1 + 1 / (1 + k)**2]
        #
        # Now, prod(1 + 1/(1 + k)**2, k=0..inf) = sinh(pi)/(2*pi) < inf.
        def eta_strategy(k, x, F):
            return f_k / (1 + k)**2

    func, x0, x0_shape, F_k, is_complex = _wrap_func(func, x0, args)

    def f(x):
        if nfev[0] >= maxfev:
            raise _NoConvergence()

        F = func(x)
        f = np.linalg.norm(F)**nexp
        nfev[0] += 1
        return f, F

    nfev = [1]

    k = 0
    x_k = x0
    f_k = np.linalg.norm(F_k)**nexp
    sigma_k = 1.0

    if fnorm is None:
        def fnorm(F):
            return f_k**(1.0/nexp)

    prev_fs = [f_k]
    F_0_norm = fnorm(F_k)

    converged = False
    message = "too many function evaluations required"

    while True:
        F_k_norm = fnorm(F_k)

        if disp:
            print("iter %d: ||F|| = %g, sigma = %g" % (k, F_k_norm, sigma_k))

        if callback is not None:
            callback(x_k, F_k)

        if F_k_norm < ftol * F_0_norm + fatol:
            # Converged!
            message = "successful convergence"
            converged = True
            break

        # Control spectral parameter, from [2]
        if abs(sigma_k) > 1/sigma_eps:
            sigma_k = 1/sigma_eps * np.sign(sigma_k)
        elif abs(sigma_k) < sigma_eps:
            sigma_k = sigma_eps

        # Line search direction
        d = -sigma_k * F_k

        # Nonmonotone line search
        eta = eta_strategy(k, x_k, F_k)
        try:
            alpha, xp, fp, Fp = _nonmonotone_line_search_cruz(f, x_k, d, prev_fs, eta=eta)
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
        prev_fs.append(fp)
        if len(prev_fs) > M:
            prev_fs.pop(0)

        k += 1

    x = _wrap_complex_result(x_k, is_complex, shape=x0_shape)
    F = _wrap_complex_result(F_k, is_complex)

    result = OptimizeResult(x=x, success=converged,
                            message=message,
                            fun=F, nfev=nfev[0], nit=k)

    return result


def _real2complex(x):
    return np.ascontiguousarray(x, dtype=float).view(np.complex128)


def _complex2real(z):
    return np.ascontiguousarray(z, dtype=complex).view(np.float64)


def _wrap_func(func, x0, args=()):
    x0 = np.asarray(x0)
    x0_shape = x0.shape
    F = np.asarray(func(x0, *args)).ravel()
    is_complex = np.iscomplexobj(x0) or np.iscomplexobj(F)
    x0 = x0.ravel()

    if is_complex:
        def wrap_func(x):
            z = _real2complex(x).reshape(x0_shape)
            v = np.asarray(func(z, *args)).ravel()
            return _complex2real(v)

        x0 = _complex2real(x0)
        F = _complex2real(F)
    else:
        def wrap_func(x, *a, **kw):
            x = x.reshape(x0_shape)
            return np.asarray(func(x, *args)).ravel()

    return wrap_func, x0, x0_shape, F, is_complex


def _wrap_complex_result(result, is_complex, shape=None):
    if is_complex:
        z = _real2complex(result)
    else:
        z = result
    if shape is not None:
        z = z.reshape(shape)
    return z
