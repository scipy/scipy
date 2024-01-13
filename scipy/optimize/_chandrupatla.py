import numpy as np
from scipy.optimize._optimize import OptimizeResult
from scipy.optimize._zeros_py import (_scalar_optimization_initialize,  # noqa: F401
                                      _chandrupatla_iv,
                                      _scalar_optimization_loop,
                                      _ECONVERGED, _ESIGNERR, _ECONVERR,
                                      _EVALUEERR, _ECALLBACK, _EINPROGRESS)


def _chandrupatla_minimize(func, x1, x2, x3, *, args=(), xatol=None,
                           xrtol=None, fatol=None, frtol=None, maxiter=100,
                           callback=None):
    """Find the minimizer of an elementwise function.

    For each element of the output of `func`, `_chandrupatla_minimize` seeks
    the scalar minimizer that minimizes the element. This function allows for
    `x1`, `x2`, `x3`, and the elements of `args` to be arrays of any
    broadcastable shapes.

    Parameters
    ----------
    func : callable
        The function whose minimizer is desired. The signature must be::

            func(x: ndarray, *args) -> ndarray

         where each element of ``x`` is a finite real and ``args`` is a tuple,
         which may contain an arbitrary number of arrays that are broadcastable
         with `x`. ``func`` must be an elementwise function: each element
         ``func(x)[i]`` must equal ``func(x[i])`` for all indices ``i``.
         `_chandrupatla` seeks an array ``x`` such that ``func(x)`` is an array
         of minima.
    x1, x2, x3 : array_like
        The abscissae of a standard scalar minimization bracket. A bracket is
        valid if ``x1 < x2 < x3`` and ``func(x1) > func(x2) <= func(x3)``.
        Must be broadcastable with one another and `args`.
    args : tuple, optional
        Additional positional arguments to be passed to `func`.  Must be arrays
        broadcastable with `x1`, `x2`, and `x3`. If the callable to be
        differentiated requires arguments that are not broadcastable with `x`,
        wrap that callable with `func` such that `func` accepts only `x` and
        broadcastable arrays.
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
        attributes. (The descriptions are written as though the values will be
        scalars; however, if `func` returns an array, the outputs will be
        arrays of the same shape.)

        success : bool
            ``True`` when the algorithm terminated successfully (status ``0``).
        status : int
            An integer representing the exit status of the algorithm.
            ``0`` : The algorithm converged to the specified tolerances.
            ``-1`` : The algorithm encountered an invalid bracket.
            ``-2`` : The maximum number of iterations was reached.
            ``-3`` : A non-finite value was encountered.
            ``-4`` : Iteration was terminated by `callback`.
            ``1`` : The algorithm is proceeding normally (in `callback` only).
        x : float
            The minimizer of the function, if the algorithm terminated
            successfully.
        fun : float
            The value of `func` evaluated at `x`.
        nfev : int
            The number of points at which `func` was evaluated.
        nit : int
            The number of iterations of the algorithm that were performed.
        xl, xm, xr : float
            The final three-point bracket.
        fl, fm, fr : float
            The function value at the bracket points.

    Notes
    -----
    Implemented based on Chandrupatla's original paper [1]_.

    If ``x1 < x2 < x3`` are the points of the bracket and ``f1 > f2 <= f3``
    are the values of ``func`` at those points, then the algorithm is
    considered to have converged when ``x3 - x1 <= abs(x2)*xrtol + xatol``
    or ``(f1 - 2*f2 + f3)/2 <= abs(f2)*frtol + fatol``. Note that first of
    these differs from the termination conditions described in [1]_. The
    default values of `xrtol` is the square root of the precision of the
    appropriate dtype, and ``xatol=fatol = frtol`` is the smallest normal
    number of the appropriate dtype.

    References
    ----------
    .. [1] Chandrupatla, Tirupathi R. (1998).
        "An efficient quadratic fit-sectioning algorithm for minimization
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
    res = _chandrupatla_iv(func, args, xatol, xrtol,
                           fatol, frtol, maxiter, callback)
    func, args, xatol, xrtol, fatol, frtol, maxiter, callback = res

    # Initialization
    xs = (x1, x2, x3)
    temp = _scalar_optimization_initialize(func, xs, args)
    xs, fs, args, shape, dtype = temp  # line split for PEP8
    x1, x2, x3 = xs
    f1, f2, f3 = fs
    phi = dtype.type(0.5 + 0.5*5**0.5)  # golden ratio
    status = np.full_like(x1, _EINPROGRESS, dtype=int)  # in progress
    nit, nfev = 0, 3  # three function evaluations performed above
    fatol = np.finfo(dtype).tiny if fatol is None else fatol
    frtol = np.finfo(dtype).tiny if frtol is None else frtol
    xatol = np.finfo(dtype).tiny if xatol is None else xatol
    xrtol = np.sqrt(np.finfo(dtype).eps) if xrtol is None else xrtol

    # Ensure that x1 < x2 < x3 initially.
    xs, fs = np.vstack((x1, x2, x3)), np.vstack((f1, f2, f3))
    i = np.argsort(xs, axis=0)
    x1, x2, x3 = np.take_along_axis(xs, i, axis=0)
    f1, f2, f3 = np.take_along_axis(fs, i, axis=0)
    q0 = x3.copy()  # "At the start, q0 is set at x3..." ([1] after (7))

    work = OptimizeResult(x1=x1, f1=f1, x2=x2, f2=f2, x3=x3, f3=f3, phi=phi,
                          xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol,
                          nit=nit, nfev=nfev, status=status, q0=q0, args=args)
    res_work_pairs = [('status', 'status'),
                      ('x', 'x2'), ('fun', 'f2'),
                      ('nit', 'nit'), ('nfev', 'nfev'),
                      ('xl', 'x1'), ('xm', 'x2'), ('xr', 'x3'),
                      ('fl', 'f1'), ('fm', 'f2'), ('fr', 'f3')]

    def pre_func_eval(work):
        # `_check_termination` is called first -> `x3 - x2 > x2 - x1`
        # But let's calculate a few terms that we'll reuse
        x21 = work.x2 - work.x1
        x32 = work.x3 - work.x2

        # [1] Section 3. "The quadratic minimum point Q1 is calculated using
        # the relations developed in the previous section." [1] Section 2 (5/6)
        A = x21 * (work.f3 - work.f2)
        B = x32 * (work.f1 - work.f2)
        C = A / (A + B)
        # q1 = C * (work.x1 + work.x2) / 2 + (1 - C) * (work.x2 + work.x3) / 2
        q1 = 0.5 * (C*(work.x1 - work.x3) + work.x2 + work.x3)  # much faster
        # this is an array, so multiplying by 0.5 does not change dtype

        # "If Q1 and Q0 are sufficiently close... Q1 is accepted if it is
        # sufficiently away from the inside point x2"
        i = abs(q1 - work.q0) < 0.5 * abs(x21)  # [1] (7)
        xi = q1[i]
        # Later, after (9), "If the point Q1 is in a +/- xtol neighborhood of
        # x2, the new point is chosen in the larger interval at a distance
        # tol away from x2."
        # See also QBASIC code after "Accept Ql adjust if close to X2".
        j = abs(q1[i] - work.x2[i]) <= work.xtol[i]
        xi[j] = work.x2[i][j] + np.sign(x32[i][j]) * work.xtol[i][j]

        # "If condition (7) is not satisfied, golden sectioning of the larger
        # interval is carried out to introduce the new point."
        # (For simplicity, we go ahead and calculate it for all points, but we
        # change the elements for which the condition was satisfied.)
        x = work.x2 + (2 - work.phi) * x32
        x[i] = xi

        # "We define Q0 as the value of Q1 at the previous iteration."
        work.q0 = q1
        return x

    def post_func_eval(x, f, work):
        # Standard logic for updating a three-point bracket based on a new
        # point. In QBASIC code, see "IF SGN(X-X2) = SGN(X3-X2) THEN...".
        # There is an awful lot of data copying going on here; this would
        # probably benefit from code optimization or implementation in Pythran.
        i = np.sign(x - work.x2) == np.sign(work.x3 - work.x2)
        xi, x1i, x2i, x3i = x[i], work.x1[i], work.x2[i], work.x3[i],
        fi, f1i, f2i, f3i = f[i], work.f1[i], work.f2[i], work.f3[i]
        j = fi > f2i
        x3i[j], f3i[j] = xi[j], fi[j]
        j = ~j
        x1i[j], f1i[j], x2i[j], f2i[j] = x2i[j], f2i[j], xi[j], fi[j]

        ni = ~i
        xni, x1ni, x2ni, x3ni = x[ni], work.x1[ni], work.x2[ni], work.x3[ni],
        fni, f1ni, f2ni, f3ni = f[ni], work.f1[ni], work.f2[ni], work.f3[ni]
        j = fni > f2ni
        x1ni[j], f1ni[j] = xni[j], fni[j]
        j = ~j
        x3ni[j], f3ni[j], x2ni[j], f2ni[j] = x2ni[j], f2ni[j], xni[j], fni[j]

        work.x1[i], work.x2[i], work.x3[i] = x1i, x2i, x3i
        work.f1[i], work.f2[i], work.f3[i] = f1i, f2i, f3i
        work.x1[ni], work.x2[ni], work.x3[ni] = x1ni, x2ni, x3ni,
        work.f1[ni], work.f2[ni], work.f3[ni] = f1ni, f2ni, f3ni

    def check_termination(work):
        # Check for all terminal conditions and record statuses.
        stop = np.zeros_like(work.x1, dtype=bool)  # termination condition met

        # Bracket is invalid; stop and don't return minimizer/minimum
        i = ((work.f2 > work.f1) | (work.f2 > work.f3))
        work.x2[i], work.f2[i] = np.nan, np.nan
        stop[i], work.status[i] = True, _ESIGNERR

        # Non-finite values; stop and don't return minimizer/minimum
        finite = np.isfinite(work.x1+work.x2+work.x3+work.f1+work.f2+work.f3)
        i = ~(finite | stop)
        work.x2[i], work.f2[i] = np.nan, np.nan
        stop[i], work.status[i] = True, _EVALUEERR

        # [1] Section 3 "Points 1 and 3 are interchanged if necessary to make
        # the (x2, x3) the larger interval."
        # Note: I had used np.choose; this is much faster. This would be a good
        # place to save e.g. `work.x3 - work.x2` for reuse, but I tried and
        # didn't notice a speed boost, so let's keep it simple.
        i = abs(work.x3 - work.x2) < abs(work.x2 - work.x1)
        temp = work.x1[i]
        work.x1[i] = work.x3[i]
        work.x3[i] = temp
        temp = work.f1[i]
        work.f1[i] = work.f3[i]
        work.f3[i] = temp

        # [1] Section 3 (bottom of page 212)
        # "We set a tolerance value xtol..."
        work.xtol = abs(work.x2) * work.xrtol + work.xatol  # [1] (8)
        # "The convergence based on interval is achieved when..."
        # Note: Equality allowed in case of `xtol=0`
        i = abs(work.x3 - work.x2) <= 2 * work.xtol  # [1] (9)

        # "We define ftol using..."
        ftol = abs(work.f2) * work.frtol + work.fatol  # [1] (10)
        # "The convergence based on function values is achieved when..."
        # Note 1: modify in place to incorporate tolerance on function value.
        # Note 2: factor of 2 is not in the text; see QBASIC start of DO loop
        i |= (work.f1 - 2 * work.f2 + work.f3) <= 2*ftol  # [1] (11)
        i &= ~stop
        stop[i], work.status[i] = True, _ECONVERGED

        return stop

    def post_termination_check(work):
        pass

    def customize_result(res, shape):
        xl, xr, fl, fr = res['xl'], res['xr'], res['fl'], res['fr']
        i = res['xl'] < res['xr']
        res['xl'] = np.choose(i, (xr, xl))
        res['xr'] = np.choose(i, (xl, xr))
        res['fl'] = np.choose(i, (fr, fl))
        res['fr'] = np.choose(i, (fl, fr))
        return shape

    return _scalar_optimization_loop(work, callback, shape,
                                     maxiter, func, args, dtype,
                                     pre_func_eval, post_func_eval,
                                     check_termination, post_termination_check,
                                     customize_result, res_work_pairs)
