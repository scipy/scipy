import numpy as np
from ._optimize import OptimizeResult, _call_callback_maybe_halt


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


def _chandrupatla_minimize(func, a, b, *, args=(), xatol=_xtol, xrtol=_rtol,
                           fatol=None, frtol=0, maxiter=_iter, callback=None):
    """Find the root of an elementwise function using Chandrupatla's algorithm.

    For each element of the output of `func`, `chandrupatla` seeks the scalar
    root that makes the element 0. This function allows for `a`, `b`, and the
    output of `func` to be of any broadcastable shapes.

    Parameters
    ----------
    func : callable
        The function whose root is desired. The signature must be::

            func(x: ndarray, *args) -> ndarray

         where each element of ``x`` is a finite real and ``args`` is a tuple,
         which may contain an arbitrary number of components of any type(s).
         ``func`` must be an elementwise function: each element ``func(x)[i]``
         must equal ``func(x[i])`` for all indices ``i``. `_chandrupatla`
         seeks an array ``x`` such that ``func(x)`` is an array of zeros.
    a, b : array_like
        The lower and upper bounds of the root of the function. Must be
        broadcastable with one another.
    args : tuple, optional
        Additional positional arguments to be passed to `func`.
    xatol, xrtol, fatol, frtol : float, optional
        Absolute and relative tolerances on the root and function value.
        See Notes for details.
    maxiter : int, optional
        The maximum number of iterations of the algorithm to perform.
    callback : callable, optional
        An optional user-supplied function to be called before the first
        iteration and after each iteration.
        Called as ``callback(res)``, where ``res`` is an ``OptimizeResult``
        similar to that returned by `_chandrupatla` (but containing the current
        iterate's values of all variables). If `callback` raises a
        ``StopIteration``, the algorithm will terminate immediately and
        `_chandrupatla` will return a result.

    Returns
    -------
    res : OptimizeResult
        An instance of `scipy.optimize.OptimizeResult` with the following
        attributes. The descriptions are written as though the values will be
        scalars; however, if `func` returns an array, the outputs will be
        arrays of the same shape.

        x : float
            The root of the function, if the algorithm terminated successfully.
        nfev : int
            The number of times the function was called to find the root.
        nit : int
            The number of iterations of Chandrupatla's algorithm performed.
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
        fun : float
            The value of `func` evaluated at `x`.
        xl, xr : float
            The lower and upper ends of the bracket.
        fl, fr : float
            The function value at the lower and upper ends of the bracket.

    Notes
    -----
    Implemented based on Chandrupatla's original paper [1]_.

    If ``xl`` and ``xr`` are the left and right ends of the bracket,
    ``xmin = xl if abs(func(xl)) <= abs(func(xr)) else xr``,
    and ``fmin0 = min(func(a), func(b))``, then the algorithm is considered to
    have converged when ``abs(xr - xl) < xatol + abs(xmin) * xrtol`` or
    ``fun(xmin) <= fatol + abs(fmin0) * frtol``. This is equivalent to the
    termination condition described in [1]_ with ``xrtol = 4e-10``,
    ``xatol = 1e-5``, and ``fatol = frtol = 0``. The default values are
    ``xatol = 2e-12``, ``xrtol = 4 * np.finfo(float).eps``, ``frtol = 0``,
    and ``fatol`` is the smallest normal number of the ``dtype`` returned
    by ``func``.

    References
    ----------

    .. [1] Chandrupatla, Tirupathi R.
        "A new hybrid quadratic/bisection algorithm for finding the zero of a
        nonlinear function without using derivatives".
        Advances in Engineering Software, 28(3), 145-149.
        https://doi.org/10.1016/s0965-9978(96)00051-8

    See Also
    --------
    brentq, brenth, ridder, bisect, newton

    Examples
    --------
    >>> from scipy import optimize
    >>> def f(x, c):
    ...     return x**3 - 2*x - c
    >>> c = 5
    >>> res = optimize._zeros_py._chandrupatla(f, 0, 3, args=(c,))
    >>> res.x
    2.0945514818937463

    >>> c = [3, 4, 5]
    >>> res = optimize._zeros_py._chandrupatla(f, 0, 3, args=(c,))
    >>> res.x
    array([1.8932892 , 2.        , 2.09455148])

    """
    res = _chandrupatla_iv(func, a, b, args, xatol, xrtol,
                           fatol, frtol, maxiter, callback)
    func, a, b, args, xatol, xrtol, fatol, frtol, maxiter, callback = res

    # Initialization
    x1, f1, x2, f2, shape, dtype = _chandrupatla_initialize(func, a, b, args)
    status = np.full_like(x1, _EINPROGRESS, dtype=int)  # in progress
    active = np.arange(x1.size)  # indices of in-progress elements
    nit, nfev = 0, 2  # two function evaluations performed above
    cb_terminate = False
    fatol = np.finfo(dtype).tiny if fatol is None else fatol
    frtol = frtol * np.minimum(np.abs(f1), np.abs(f2))
    tols = dict(xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol)
    # Elements of `x1`, `f1`, etc., are stored in this `OptimizeResult`
    # once a termination condition is met, and then the arrays are condensed
    # to reduce unnecessary computation.
    res = OptimizeResult(x=x1.copy(), fun=f1.copy(), xl=x1.copy(),
                         fl=f1.copy(), xr=x2.copy(), fr=f2.copy(),
                         nit=np.full_like(status, nit)[()],
                         nfev=np.full_like(status, nfev)[()],
                         status=status.copy(), success=(status==0))

    temp = _chandrupatla_check_termination(x1, f1, x2, f2, None, None, res,
                                           active, status, nfev, nit, tols)
    xmin, fmin, x1, f1, x2, f2, x3, f3, active, status, tl = temp

    if callback is not None:
        temp = _chandrupatla_prepare_result(shape, res, active, nit, nfev)
        if _call_callback_maybe_halt(callback, temp):
            cb_terminate = True

    t = 0.5
    while nit < maxiter and active.size and not cb_terminate:
        # [1] Figure 1 (first box)
        x = x1 + t * (x2 - x1)

        # For now, we assume that `func` requires arguments of the original
        # shapes. However, `x` is compressed (contains only active elements).
        # Therefore, we inflate `x` by creating a full-size array, copying the
        # elements of `x` into it, and making it the expected shape.
        # TODO: allow user to specify that `func` works with compressed input
        x_full = res.x.copy()
        x_full[active] = x
        x_full = x_full.reshape(shape)
        f = func(x_full, *args)
        nfev += 1
        # Ensure that the outpuf of `func` is an array of the appropriate
        # dtype, ravel it (because we work with 1D arrays for simplicity), and
        # compress it to contain only active elements.
        f = np.asarray(f, dtype=dtype).ravel()[active]

        # [1] Figure 1 (first diamond and boxes)
        # Note: y/n are reversed in figure; compare to BASIC in appendix
        x3, f3 = x2.copy(), f2.copy()
        j = np.sign(f) == np.sign(f1)
        nj = ~j
        x3[j], f3[j] = x1[j], f1[j]
        x2[nj], f2[nj] = x1[nj], f1[nj]
        x1, f1 = x, f

        # [1] Figure 1 (second diamond)
        nit += 1
        temp = _chandrupatla_check_termination(x1, f1, x2, f2, x3, f3, res,
                                               active, status, nfev, nit, tols)
        xmin, fmin, x1, f1, x2, f2, x3, f3, active, status, tl = temp

        if callback is not None:
            temp = _chandrupatla_prepare_result(shape, res, active, nit, nfev)
            if _call_callback_maybe_halt(callback, temp):
                cb_terminate = True
                break
        if active.size==0:
            break

        # [1] Figure 1 (third diamond and boxes / Equation 1)
        xi1 = (x1 - x2) / (x3 - x2)
        phi1 = (f1 - f2) / (f3 - f2)
        alpha = (x3 - x1) / (x2 - x1)
        j = ((1 - np.sqrt(1 - xi1)) < phi1) & (phi1 < np.sqrt(xi1))

        f1j, f2j, f3j, alphaj = f1[j], f2[j], f3[j], alpha[j]
        t = np.full_like(alpha, 0.5)
        t[j] = (f1j / (f1j - f2j) * f3j / (f3j - f2j)
                - alphaj * f1j / (f3j - f1j) * f2j / (f2j - f3j))

        # [1] Figure 1 (last box; see also BASIC in appendix with comment
        # "Adjust T Away from the Interval Boundary")
        t = np.clip(t, tl, 1-tl)

    res.status[active] = _ECALLBACK if cb_terminate else _ECONVERR
    return _chandrupatla_prepare_result(shape, res, active, nit, nfev)

def _chandrupatla_iv(func, a, b, args, xatol, xrtol,
                     fatol, frtol, maxiter, callback):

    if not callable(func):
        raise ValueError('`func` must be callable.')

    # a and b have more complex IV that is taken care of during initialization

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

    return func, a, b, args, xatol, xrtol, fatol, frtol, maxiter, callback

def _chandrupatla_initialize(func, a, b, args):
    # initializing left and right bracket and function value arrays

    # Try to preserve `dtype`, but we need to ensure that the arguments are at
    # least floats before passing them into the function because integers
    # can overflow and cause failure.
    x1, x2 = np.broadcast_arrays(a, b)  # broadcast and rename
    xt = np.result_type(x1.dtype, x2.dtype)
    xt = np.float64 if np.issubdtype(xt, np.integer) else xt
    x1, x2 = x1.astype(xt, copy=False)[()], x2.astype(xt, copy=False)[()]
    f1, f2 = func(x1, *args), func(x2, *args)

    # It's possible that the functions will return multiple outputs for each
    # scalar input, so we need to broadcast again. All arrays need to be,
    # writable, so we'll need to copy after broadcasting. Finally, we're going
    # to be doing operations involving `x1`, `x2`, `f1`, and `f2` throughout,
    # so figure out the right type from the outset.
    x1, f1, x2, f2 = np.broadcast_arrays(x1, f1, x2, f2)
    ft = np.result_type(x1.dtype, x2.dtype, f1.dtype, f2.dtype)
    ft = np.float64 if np.issubdtype(ft, np.integer) else ft
    if not np.issubdtype(np.result_type(x1, x2, f1, f2), np.floating):
        raise ValueError("Bracket and function output must be real numbers.")
    x1, f1, x2, f2 = x1.astype(ft), f1.astype(ft), x2.astype(ft), f2.astype(ft)

    # To ensure that we can do indexing, we'll work with at least 1d arrays,
    # but remember the appropriate shape of the output.
    shape = x1.shape
    x1, f1, x2, f2, = x1.ravel(), f1.ravel(), x2.ravel(), f2.ravel()
    return x1, f1, x2, f2, shape, ft


def _chandrupatla_check_termination(x1, f1, x2, f2, x3, f3, res,
                                    active, status, nfev, nit, tols):
    # Check for all terminal conditions and record statuses.

    # See [1] Section 4 (first two sentences)
    i = np.abs(f1) < np.abs(f2)
    xmin = np.choose(i, (x2, x1))
    fmin = np.choose(i, (f2, f1))
    stop = np.zeros_like(x1, dtype=bool)  # termination condition met

    # This is the convergence criterion used in bisect. Chandrupatla's
    # criterion is equivalent to this except with a factor of 4 on `xrtol`.
    dx = abs(x2 - x1)
    tol = abs(xmin)*tols['xrtol'] + tols['xatol']
    i = dx < tol
    # Modify in place to incorporate tolerance on function value. Note that
    # `frtol` has been redefined as `frtol = frtol * np.minimum(f1, f2)`,
    # where `f1` and `f2` are the function evaluated at the original ends of
    # the bracket.
    i |= np.abs(fmin) <= tols['fatol'] + tols['frtol']
    stop[i], status[i] = True, _ECONVERGED

    i = (np.sign(f1) == np.sign(f2)) & ~stop
    xmin[i], fmin[i], stop[i], status[i] = np.nan, np.nan, True, _ESIGNERR

    i = ~((np.isfinite(x1) & np.isfinite(x2)
           & np.isfinite(f1) & np.isfinite(f2)) | stop)
    xmin[i], fmin[i], stop[i], status[i] = np.nan, np.nan, True, _EVALUEERR

    if np.any(stop):
        # update the result object with the elements for which termination
        # condition has been met
        active_stop = active[stop]
        res.x[active_stop] = xmin[stop]
        res.fun[active_stop] = fmin[stop]
        res.xl[active_stop] = x1[stop]
        res.xr[active_stop] = x2[stop]
        res.fl[active_stop] = f1[stop]
        res.fr[active_stop] = f2[stop]
        res.status[active_stop] = status[stop]
        res.nfev[active_stop] = nfev
        res.nit[active_stop] = nit
        res.success[active_stop] = res.status[active_stop] == 0

        # compress the arrays to avoid unnecessary computation
        proceed = ~stop
        active = active[proceed]
        xmin = xmin[proceed]
        fmin = fmin[proceed]
        x1 = x1[proceed]
        f1 = f1[proceed]
        x2 = x2[proceed]
        f2 = f2[proceed]
        x3 = x3[proceed] if x3 is not None else x3
        f3 = f3[proceed] if f3 is not None else f3
        status = status[proceed]
        tols['frtol'] = tols['frtol'][proceed]
        tol = tol[proceed]
        dx = dx[proceed]

    # See [1] Appendix, "Adjust T Away from the Interval Boundary"
    tl = 0.5 * tol / dx

    return xmin, fmin, x1, f1, x2, f2, x3, f3, active, status, tl


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
