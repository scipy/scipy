import numpy as np
from scipy.optimize._optimize import OptimizeResult, _call_callback_maybe_halt
from scipy.optimize._zeros_py import _scalar_optimization_initialize

_iter = 100
_xtol = 2e-12
_rtol = 4 * np.finfo(float).eps

__all__ = []

# Must agree with CONVERGED, SIGNERR, CONVERR, ...  in zeros.h
_ECONVERGED = 0
_ESIGNERR = -1
_ECONVERR = -2
_EVALUEERR = -3
_ECALLBACK = -4
_EINPROGRESS = 1

# TODO:
#  - Figure out whether we want to follow original termination conditions
#  - In this and `_chandrupatla`, fix result object not updating
#  - Check over object returned by `_chandrupatla_minimize`. Make sure that
#    attributes are printed in a reasonable order and that `xl < x < xr`.
#    Should `xm` even be reported if it's identical to `x`?

def _chandrupatla_minimize(func, x1, x2, x3, *, args=(), xatol=_xtol,
                           xrtol=_rtol, fatol=None, frtol=0, maxiter=_iter,
                           callback=None):
    """Find the minimizer of an elementwise function.

    For each element of the output of `func`, `_chandrupatla_minimize` seeks
    the scalar minimizer that minimizes the element. This function allows for
    `x1`, `x2`, `x3`, and the output of `func` to be of any broadcastable
    shapes.

    Parameters
    ----------
    func : callable
        The function whose minimizer is desired. The signature must be::

            func(x: ndarray, *args) -> ndarray

         where each element of ``x`` is a finite real and ``args`` is a tuple,
         which may contain an arbitrary number of components of any type(s).
         ``func`` must be an elementwise function: each element ``func(x)[i]``
         must equal ``func(x[i])`` for all indices ``i``. `_chandrupatla`
         seeks an array ``x`` such that ``func(x)`` is an array of minima.
    x1, x2, x3 : array_like
        The abscissae of a standard scalar minimization bracket. A bracket is
        valid if ``x1 < x2 < x3`` and ``func(x1) > func(x2) < func(x3)``.
        Must be broadcastable with one another.
    args : tuple, optional
        Additional positional arguments to be passed to `func`.
    xatol, xrtol, fatol, frtol : float, optional
        Absolute and relative tolerances on the minimizer and function value.
        See Notes for details.
    maxiter : int, optional
        The maximum number of iterations of the algorithm to perform.
    callback : callable, optional
        An optional user-supplied function to be called before the first
        iteration and after each iteration.
        Called as ``callback(res)``, where ``res`` is an ``OptimizeResult``
        similar to that returned by `_chandrupatla_minimize` (but containing
        the current iterate's values of all variables). If `callback` raises a
        ``StopIteration``, the algorithm will terminate immediately and
        `_chandrupatla_minimize` will return a result.

    Returns
    -------
    res : OptimizeResult
        An instance of `scipy.optimize.OptimizeResult` with the following
        attributes. The descriptions are written as though the values will be
        scalars; however, if `func` returns an array, the outputs will be
        arrays of the same shape.

        x : float
            The minimizer of the function, if the algorithm terminated
            successfully.
        fun : float
            The value of `func` evaluated at `x`.
        nfev : int
            The number of times the function was called.
        nit : int
            The number of iterations of the algorithm that were performed.
        status : int
            An integer representing the exit status of the algorithm.
            ``0`` : The algorithm converged to the specified tolerances.
            ``-1`` : The algorithm encountered an invalid bracket.
            ``-2`` : The maximum number of iterations was reached.
            ``-3`` : A non-finite value was encountered.
            ``-4`` : Iteration was terminated by `callback`.
            ``1`` : The algorithm is proceeding normally (in `callback` only).
        success : bool
            ``True`` when the algorithm terminated successfully (status ``0``).
        x1, x2, x3 : float
            The final three-point bracket.
        f1, f2, f3 : float
            The function value at the bracket points.

    Notes
    -----
    Implemented based on Chandrupatla's original paper [1]_.

    If ``x1 < x2 < x3`` are the points of the bracket and ``f1 > f2 <= f3``
    are the values of ``func`` at those points, then the algorithm is
    considered to have converged when ``x3 - x1 <= abs(x2)*xrtol + xatol``
    or ``(f1 - 2*f2 + f3)/2 <= abs(f2)*frtol + fatol``. Note that first of
    these differs from the termination conditions described in [1]_. The
    default values of the tolerances are ``xatol = 1e-12``,
    ``xrtol = 4 * np.finfo(float).eps``, and ``fatol = frtol = 0``.

    References
    ----------
    .. [1] Chandrupatla, Tirupathi R. (1998).
        "An efficient quadratic fit—sectioning algorithm for minimization
        without derivatives".
        Computer Methods in Applied Mechanics and Engineering, 152 (1-2),
        211-217. https://doi.org/10.1016/S0045-7825(97)00190-4

    See Also
    --------
    golden, brent, bounded

    Examples
    --------
    >>> from scipy.optimize._chandrupatla import _chandrupatla_minimize
    >>> def f(x, args=1):
    ...     return (x - args)**2
    >>> res = _chandrupatla_minimize(f, -5, 0, 5)
    >>> res.x
    1.0
    >>> c = [1, 1.5, 2]
    >>> res = _chandrupatla_minimize(f, -5, 0, 5, args=(c,))
    >>> res.x
    array([1. , 1.5, 2. ])
    """
    # _chandrupatla_iv can be used essentially verbatim. Validation of `a`/`b`
    # (now `x1`, `x2`, `x3`) happen separately
    res = _chandrupatla_iv(func, x1, x2, x3, args, xatol, xrtol,
                           fatol, frtol, maxiter, callback)
    func, x1, x2, x3, args, xatol, xrtol, fatol, frtol, maxiter, callback = res

    # Initialization
    # _chandrupatla_initialize can be used by both if generalized to accept
    # an arbitrary number of `x` elements
    xs = (x1, x2, x3)
    xs, fs, args, shape, dtype = _scalar_optimization_initialize(func, xs, args)
    x1, x2, x3 = xs
    f1, f2, f3 = fs
    q0 = x3  # At the start, q0 is set at x3...
    phi = 1.61803398875

    # Ensure that x1 < x2 < x3 initially.
    xs, fs = np.vstack((x1, x2, x3)), np.vstack((f1, f2, f3))
    i = np.argsort(xs, axis=0)
    x1, x2, x3 = np.take_along_axis(xs, i, axis=0)
    f1, f2, f3 = np.take_along_axis(fs, i, axis=0)

    # all this is the same except nfev should be initialized to the number
    # of different `x` arrays
    status = np.full_like(x1, _EINPROGRESS, dtype=int)  # in progress
    active = np.arange(x1.size)  # indices of in-progress elements
    nit, nfev = 0, 3  # three function evaluations performed above
    cb_terminate = False
    fatol = np.finfo(dtype).tiny if fatol is None else fatol
    tols = dict(xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol)

    # Elements of `x1`, `f1`, etc., are stored in this `OptimizeResult`
    # once a termination condition is met, and then the arrays are condensed
    # to reduce unnecessary computation.
    res = OptimizeResult(x=x1.copy(), fun=f1.copy(), xl=x1.copy(),
                         fl=f1.copy(), xm=x2.copy(), fm=f2.copy(),
                         xr=x3.copy(), fr=f3.copy(),
                         nit=np.full_like(status, nit)[()],
                         nfev=np.full_like(status, nfev)[()],
                         status=status.copy(), success=(status==0))

    temp = _chandrupatla_check_termination(x1, f1, x2, f2, x3, f3, q0, res,
                                           active, status, nfev, nit, tols)
    x1, f1, x2, f2, x3, f3, q0, xtol, active, status = temp

    if callback is not None:
        temp = _chandrupatla_prepare_result(shape, res, active, nit, nfev)
        if _call_callback_maybe_halt(callback, temp):
            cb_terminate = True

    while nit < maxiter and active.size and not cb_terminate:

        A = (x2 - x1) * (f3 - f2)
        B = (x3 - x2) * (f1 - f2)
        C = A / (A + B)
        q1 = C * (x1 + x2) / 2 + (1 - C) * (x2 + x3) / 2

        dq = abs(q1 - q0)
        dx12 = abs(x2 - x1)

        i = dq < 0.5 * dx12
        x = x2 + (2 - phi) * (x3 - x2)

        j = abs(q1[i] - x2[i]) > xtol[i]
        xi = x[i]
        xi[j] = q1[i][j]
        xi[~j] = x2[i][~j] + np.sign(x3[i][~j] - x2[i][~j]) * xtol[i][~j]
        x[i] = xi

        q0 = q1
        x_full = res.x.copy()
        x_full[active] = x
        x_full = x_full.reshape(shape)
        args_full = [arg.reshape(shape) for arg in args]
        f = func(x_full, *args_full)
        nfev += 1
        f = np.asarray(f, dtype=dtype).ravel()[active]

        i = np.sign(x - x2) == np.sign(x3 - x2)

        # TODO: tame this mess
        xi, fi, x1i, f1i, x2i, f2i, x3i, f3i = x[i], f[i], x1[i], f1[i], x2[i], f2[i], x3[i], f3[i]
        j = fi > f2i
        x3i[j], f3i[j] = xi[j], fi[j]
        x1i[~j], f1i[~j], x2i[~j], f2i[~j] = x2i[~j], f2i[~j], xi[~j], fi[~j]

        xni, fni, x1ni, f1ni, x2ni, f2ni, x3ni, f3ni = x[~i], f[~i], x1[~i], f1[~i], x2[~i], f2[~i], x3[~i], f3[~i]
        j = fni > f2ni
        x1ni[j], f1ni[j] = xni[j], fni[j]
        x3ni[~j], f3ni[~j], x2ni[~j], f2ni[~j] = x2ni[~j], f2ni[~j], xni[~j], fni[~j]

        x[i], f[i], x1[i], f1[i], x2[i], f2[i], x3[i], f3[i] = xi, fi, x1i, f1i, x2i, f2i, x3i, f3i
        x[~i], f[~i], x1[~i], f1[~i], x2[~i], f2[~i], x3[~i], f3[~i] = xni, fni, x1ni, f1ni, x2ni, f2ni, x3ni, f3ni

        # [1] Figure 1 (second diamond)
        nit += 1
        temp = _chandrupatla_check_termination(x1, f1, x2, f2, x3, f3, q0, res,
                                               active, status, nfev, nit, tols)
        x1, f1, x2, f2, x3, f3, q0, xtol, active, status = temp

        if callback is not None:
            temp = _chandrupatla_prepare_result(shape, res, active, nit, nfev)
            if _call_callback_maybe_halt(callback, temp):
                cb_terminate = True
                break
        if active.size==0:
            break

    res.status[active] = _ECALLBACK if cb_terminate else _ECONVERR
    return _chandrupatla_prepare_result(shape, res, active, nit, nfev)


def _chandrupatla_iv(func, x1, x2, x3, args, xatol, xrtol,
                     fatol, frtol, maxiter, callback):

    if not callable(func):
        raise ValueError('`func` must be callable.')

    # x1, x2, x3 have more complex IV that is taken care of during
    # initialization

    if not np.iterable(args):
        args = (args,)

    tols = np.asarray([xatol, xrtol, fatol if fatol is not None else 1, frtol])
    if (not np.issubdtype(tols.dtype, np.number)
            or np.any(tols < 0)
            or tols.shape != (4,)):
        raise ValueError('Tolerances must be non-negative scalars.')

    maxiter_int = int(maxiter)
    if maxiter != maxiter_int or maxiter < 0:
        raise ValueError('`maxiter` must be a non-negative integer.')

    if callback is not None and not callable(callback):
        raise ValueError('`callback` must be callable.')

    return func, x1, x2, x3, args, xatol, xrtol, fatol, frtol, maxiter, callback


def _chandrupatla_check_termination(x1, f1, x2, f2, x3, f3, q0, res,
                                    active, status, nfev, nit, tols):
    # Check for all terminal conditions and record statuses.

    # See [1] Section 4 (first two sentences)
    i = abs(x2 - x1) < abs(x3 - x2)
    x1, x3 = np.choose(i, (x3, x1)), np.choose(i, (x1, x3))
    f1, f3 = np.choose(i, (f3, f1)), np.choose(i, (f1, f3))
    stop = np.zeros_like(x1, dtype=bool)  # termination condition met

    i = ((f2 > f1) | (f2 > f3))
    x2[i], f2[i], stop[i], status[i] = np.nan, np.nan, True, _ESIGNERR

    xtol = abs(x2)*tols['xrtol'] + tols['xatol']
    i = abs(x3 - x2)/2 < xtol
    # Modify in place to incorporate tolerance on function value.
    i |= (f1 - 2 * f2 + f3)/2 < abs(f2)*tols['frtol'] + tols['fatol']
    i &=  ~stop
    stop[i], status[i] = True, _ECONVERGED

    i = ~((np.isfinite(x1) & np.isfinite(x2) & np.isfinite(x3)
           & np.isfinite(f1) & np.isfinite(f2) & np.isfinite(f3)) | stop)
    x2[i], x2[i], stop[i], status[i] = np.nan, np.nan, True, _EVALUEERR

    ### This stuff can be put into a function and reused
    if np.any(stop):
        # update the result object with the elements for which termination
        # condition has been met
        active_stop = active[stop]
        res.x[active_stop] = x2[stop]
        res.fun[active_stop] = f2[stop]
        res.xl[active_stop] = x1[stop]
        res.xm[active_stop] = x2[stop]
        res.xr[active_stop] = x3[stop]
        res.fl[active_stop] = f1[stop]
        res.fm[active_stop] = f2[stop]
        res.fr[active_stop] = f3[stop]
        res.status[active_stop] = status[stop]
        res.nfev[active_stop] = nfev
        res.nit[active_stop] = nit
        res.success[active_stop] = res.status[active_stop] == 0

        # compress the arrays to avoid unnecessary computation
        proceed = ~stop
        active = active[proceed]
        x1 = x1[proceed]
        f1 = f1[proceed]
        x2 = x2[proceed]
        f2 = f2[proceed]
        x3 = x3[proceed]
        f3 = f3[proceed]
        q0 = q0[proceed]
        xtol = xtol[proceed]
        status = status[proceed]

    return x1, f1, x2, f2, x3, f3, q0, xtol, active, status


def _chandrupatla_prepare_result(shape, res, active, nit, nfev):
    res = res.copy()
    xl, xr, fl, fr = res['xl'], res['xr'], res['fl'], res['fr']
    i = res['xl'] < res['xr']
    res['xl'] = np.choose(i, (xr, xl))
    res['xr'] = np.choose(i, (xl, xr))
    res['fl'] = np.choose(i, (fr, fl))
    res['fr'] = np.choose(i, (fl, fr))
    res['nit'][active] = nit
    res['nfev'][active] = nfev
    for key, val in res.items():
        res[key] = np.reshape(val, shape)[()]
    return OptimizeResult(**res)
