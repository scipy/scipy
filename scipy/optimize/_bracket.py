import numpy as np
import scipy._lib._elementwise_iterative_method as eim
from scipy._lib._util import _RichResult

_ELIMITS = -1  # used in _bracket_root
_ESTOPONESIDE = 2  # used in _bracket_root

def _bracket_root_iv(func, a, b, min, max, factor, args, maxiter):

    if not callable(func):
        raise ValueError('`func` must be callable.')

    if not np.iterable(args):
        args = (args,)

    a = np.asarray(a)[()]
    if not np.issubdtype(a.dtype, np.number) or np.iscomplex(a).any():
        raise ValueError('`a` must be numeric and real.')

    b = a + 1 if b is None else b
    min = -np.inf if min is None else min
    max = np.inf if max is None else max
    factor = 2. if factor is None else factor
    a, b, min, max, factor = np.broadcast_arrays(a, b, min, max, factor)

    if not np.issubdtype(b.dtype, np.number) or np.iscomplex(b).any():
        raise ValueError('`b` must be numeric and real.')

    if not np.issubdtype(min.dtype, np.number) or np.iscomplex(min).any():
        raise ValueError('`min` must be numeric and real.')

    if not np.issubdtype(max.dtype, np.number) or np.iscomplex(max).any():
        raise ValueError('`max` must be numeric and real.')

    if not np.issubdtype(factor.dtype, np.number) or np.iscomplex(factor).any():
        raise ValueError('`factor` must be numeric and real.')
    if not np.all(factor > 1):
        raise ValueError('All elements of `factor` must be greater than 1.')

    maxiter = np.asarray(maxiter)
    message = '`maxiter` must be a non-negative integer.'
    if (not np.issubdtype(maxiter.dtype, np.number) or maxiter.shape != tuple()
            or np.iscomplex(maxiter)):
        raise ValueError(message)
    maxiter_int = int(maxiter[()])
    if not maxiter == maxiter_int or maxiter < 0:
        raise ValueError(message)

    if not np.all((min <= a) & (a < b) & (b <= max)):
        raise ValueError('`min <= a < b <= max` must be True (elementwise).')

    return func, a, b, min, max, factor, args, maxiter


def _bracket_root(func, a, b=None, *, min=None, max=None, factor=None,
                  args=(), maxiter=1000):
    """Bracket the root of a monotonic scalar function of one variable

    This function works elementwise when `a`, `b`, `min`, `max`, `factor`, and
    the elements of `args` are broadcastable arrays.

    Parameters
    ----------
    func : callable
        The function for which the root is to be bracketed.
        The signature must be::

            func(x: ndarray, *args) -> ndarray

        where each element of ``x`` is a finite real and ``args`` is a tuple,
        which may contain an arbitrary number of arrays that are broadcastable
        with `x`. ``func`` must be an elementwise function: each element
        ``func(x)[i]`` must equal ``func(x[i])`` for all indices ``i``.
    a, b : float array_like
        Starting guess of bracket, which need not contain a root. If `b` is
        not provided, ``b = a + 1``. Must be broadcastable with one another.
    min, max : float array_like, optional
        Minimum and maximum allowable endpoints of the bracket, inclusive. Must
        be broadcastable with `a` and `b`.
    factor : float array_like, default: 2
        The factor used to grow the bracket. See notes for details.
    args : tuple, optional
        Additional positional arguments to be passed to `func`.  Must be arrays
        broadcastable with `a`, `b`, `min`, and `max`. If the callable to be
        bracketed requires arguments that are not broadcastable with these
        arrays, wrap that callable with `func` such that `func` accepts
        only `x` and broadcastable arrays.
    maxiter : int, optional
        The maximum number of iterations of the algorithm to perform.

    Returns
    -------
    res : _RichResult
        An instance of `scipy._lib._util._RichResult` with the following
        attributes. The descriptions are written as though the values will be
        scalars; however, if `func` returns an array, the outputs will be
        arrays of the same shape.

        xl, xr : float
            The lower and upper ends of the bracket, if the algorithm
            terminated successfully.
        fl, fr : float
            The function value at the lower and upper ends of the bracket.
        nfev : int
            The number of function evaluations required to find the bracket.
            This is distinct from the number of times `func` is *called*
            because the function may evaluated at multiple points in a single
            call.
        nit : int
            The number of iterations of the algorithm that were performed.
        status : int
            An integer representing the exit status of the algorithm.

            - ``0`` : The algorithm produced a valid bracket.
            - ``-1`` : The bracket expanded to the allowable limits without finding a bracket.
            - ``-2`` : The maximum number of iterations was reached.
            - ``-3`` : A non-finite value was encountered.
            - ``-4`` : Iteration was terminated by `callback`.
            - ``1`` : The algorithm is proceeding normally (in `callback` only).
            - ``2`` : A bracket was found in the opposite search direction (in `callback` only).

        success : bool
            ``True`` when the algorithm terminated successfully (status ``0``).

    Notes
    -----
    This function generalizes an algorithm found in pieces throughout
    `scipy.stats`. The strategy is to iteratively grow the bracket `(l, r)`
     until ``func(l) < 0 < func(r)``. The bracket grows to the left as follows.

    - If `min` is not provided, the distance between `b` and `l` is iteratively
      increased by `factor`.
    - If `min` is provided, the distance between `min` and `l` is iteratively
      decreased by `factor`. Note that this also *increases* the bracket size.

    Growth of the bracket to the right is analogous.

    Growth of the bracket in one direction stops when the endpoint is no longer
    finite, the function value at the endpoint is no longer finite, or the
    endpoint reaches its limiting value (`min` or `max`). Iteration terminates
    when the bracket stops growing in both directions, the bracket surrounds
    the root, or a root is found (accidentally).

    If two brackets are found - that is, a bracket is found on both sides in
    the same iteration, the smaller of the two is returned.
    If roots of the function are found, both `l` and `r` are set to the
    leftmost root.

    """  # noqa: E501
    # Todo:
    # - find bracket with sign change in specified direction
    # - Add tolerance
    # - allow factor < 1?

    callback = None  # works; I just don't want to test it
    temp = _bracket_root_iv(func, a, b, min, max, factor, args, maxiter)
    func, a, b, min, max, factor, args, maxiter = temp

    xs = (a, b)
    temp = eim._initialize(func, xs, args)
    func, xs, fs, args, shape, dtype = temp  # line split for PEP8

    # The approach is to treat the left and right searches as though they were
    # (almost) totally independent one-sided bracket searches. (The interaction
    # is considered when checking for termination and preparing the result
    # object.)
    # `x` is the "moving" end of the bracket
    x = np.concatenate(xs)
    f = np.concatenate(fs)
    n = len(x) // 2

    # `x_last` is the previous location of the moving end of the bracket. If
    # the signs of `f` and `f_last` are different, `x` and `x_last` form a
    # bracket.
    x_last = np.concatenate((x[n:], x[:n]))
    f_last = np.concatenate((f[n:], f[:n]))
    # `x0` is the "fixed" end of the bracket.
    x0 = x_last
    # We don't need to retain the corresponding function value, since the
    # fixed end of the bracket is only needed to compute the new value of the
    # moving end; it is never returned.

    min = np.broadcast_to(min, shape).astype(dtype, copy=False).ravel()
    max = np.broadcast_to(max, shape).astype(dtype, copy=False).ravel()
    limit = np.concatenate((min, max))

    factor = np.broadcast_to(factor, shape).astype(dtype, copy=False).ravel()
    factor = np.concatenate((factor, factor))

    active = np.arange(2*n)
    args = [np.concatenate((arg, arg)) for arg in args]

    # This is needed due to inner workings of `eim._loop`.
    # We're abusing it a tiny bit.
    shape = shape + (2,)

    # `d` is for "distance".
    # For searches without a limit, the distance between the fixed end of the
    # bracket `x0` and the moving end `x` will grow by `factor` each iteration.
    # For searches with a limit, the distance between the `limit` and moving
    # end of the bracket `x` will shrink by `factor` each iteration.
    i = np.isinf(limit)
    ni = ~i
    d = np.zeros_like(x)
    d[i] = x[i] - x0[i]
    d[ni] = limit[ni] - x[ni]

    status = np.full_like(x, eim._EINPROGRESS, dtype=int)  # in progress
    nit, nfev = 0, 1  # one function evaluation per side performed above

    work = _RichResult(x=x, x0=x0, f=f, limit=limit, factor=factor,
                       active=active, d=d, x_last=x_last, f_last=f_last,
                       nit=nit, nfev=nfev, status=status, args=args,
                       xl=None, xr=None, fl=None, fr=None, n=n)
    res_work_pairs = [('status', 'status'), ('xl', 'xl'), ('xr', 'xr'),
                      ('nit', 'nit'), ('nfev', 'nfev'), ('fl', 'fl'),
                      ('fr', 'fr'), ('x', 'x'), ('f', 'f'),
                      ('x_last', 'x_last'), ('f_last', 'f_last')]

    def pre_func_eval(work):
        # Initialize moving end of bracket
        x = np.zeros_like(work.x)

        # Unlimited brackets grow by `factor` by increasing distance from fixed
        # end to moving end.
        i = np.isinf(work.limit)  # indices of unlimited brackets
        work.d[i] *= work.factor[i]
        x[i] = work.x0[i] + work.d[i]

        # Limited brackets grow by decreasing the distance from the limit to
        # the moving end.
        ni = ~i  # indices of limited brackets
        work.d[ni] /= work.factor[ni]
        x[ni] = work.limit[ni] - work.d[ni]

        return x

    def post_func_eval(x, f, work):
        # Keep track of the previous location of the moving end so that we can
        # return a narrower bracket. (The alternative is to remember the
        # original fixed end, but then the bracket would be wider than needed.)
        work.x_last = work.x
        work.f_last = work.f
        work.x = x
        work.f = f

    def check_termination(work):
        stop = np.zeros_like(work.x, dtype=bool)

        # Condition 1: a valid bracket (or the root itself) has been found
        sf = np.sign(work.f)
        sf_last = np.sign(work.f_last)
        i = (sf_last == -sf) | (sf_last == 0) | (sf == 0)
        work.status[i] = eim._ECONVERGED
        stop[i] = True

        # Condition 2: the other side's search found a valid bracket.
        # (If we just found a bracket with the rightward search, we can stop
        #  the leftward search, and vice-versa.)
        # To do this, we need to set the status of the other side's search;
        # this is tricky because `work.status` contains only the *active*
        # elements, so we don't immediately know the index of the element we
        # need to set - or even if it's still there. (That search may have
        # terminated already, e.g. by reaching its `limit`.)
        # To facilitate this, `work.active` contains a unit integer index of
        # each search. Index `k` (`k < n)` and `k + n` correspond with a
        # leftward and rightward search, respectively. Elements are removed
        # from `work.active` just as they are removed from `work.status`, so
        # we use `work.active` to help find the right location in
        # `work.status`.
        # Get the integer indices of the elements that can also stop
        also_stop = (work.active[i] + work.n) % (2*work.n)
        # Check whether they are still active.
        # To start, we need to find out where in `work.active` they would
        # appear if they are indeed there.
        j = np.searchsorted(work.active, also_stop)
        # If the location exceeds the length of the `work.active`, they are
        # not there.
        j = j[j < len(work.active)]
        # Check whether they are still there.
        j = j[also_stop == work.active[j]]
        # Now convert these to boolean indices to use with `work.status`.
        i = np.zeros_like(stop)
        i[j] = True  # boolean indices of elements that can also stop
        i = i & ~stop
        work.status[i] = _ESTOPONESIDE
        stop[i] = True

        # Condition 3: moving end of bracket reaches limit
        i = (work.x == work.limit) & ~stop
        work.status[i] = _ELIMITS
        stop[i] = True

        # Condition 4: non-finite value encountered
        i = ~(np.isfinite(work.x) & np.isfinite(work.f)) & ~stop
        work.status[i] = eim._EVALUEERR
        stop[i] = True

        return stop

    def post_termination_check(work):
        pass

    def customize_result(res, shape):
        n = len(res['x']) // 2

        # Because we treat the two one-sided searches as though they were
        # independent, what we keep track of in `work` and what we want to
        # return in `res` look quite different. Combine the results from the
        # two one-sided searches before reporting the results to the user.
        # - "a" refers to the leftward search (the moving end started at `a`)
        # - "b" refers to the rightward search (the moving end started at `b`)
        # - "l" refers to the left end of the bracket (closer to -oo)
        # - "r" refers to the right end of the bracket (closer to +oo)
        xal = res['x'][:n]
        xar = res['x_last'][:n]
        xbl = res['x_last'][n:]
        xbr = res['x'][n:]

        fal = res['f'][:n]
        far = res['f_last'][:n]
        fbl = res['f_last'][n:]
        fbr = res['f'][n:]

        # Initialize the brackets and corresponding function values to return
        # to the user. Brackets may not be valid (e.g. there is no root,
        # there weren't enough iterations, NaN encountered), but we still need
        # to return something. One option would be all NaNs, but what I've
        # chosen here is the left- and right-most points at which the function
        # has been evaluated. This gives the user some information about what
        # interval of the real line has been searched and shows that there is
        # no sign change between the two ends.
        xl = xal.copy()
        fl = fal.copy()
        xr = xbr.copy()
        fr = fbr.copy()

        # `status` indicates whether the bracket is valid or not. If so,
        # we want to adjust the bracket we return to be the narrowest possible
        # given the points at which we evaluated the function.
        # For example if bracket "a" is valid and smaller than bracket "b" OR
        # if bracket "a" is valid and bracket "b" is not valid, we want to
        # return bracket "a" (and vice versa).
        sa = res['status'][:n]
        sb = res['status'][n:]

        da = xar - xal
        db = xbr - xbl

        i1 = ((da <= db) & (sa == 0)) | ((sa == 0) & (sb != 0))
        i2 = ((db <= da) & (sb == 0)) | ((sb == 0) & (sa != 0))

        xr[i1] = xar[i1]
        fr[i1] = far[i1]
        xl[i2] = xbl[i2]
        fl[i2] = fbl[i2]

        # Finish assembling the result object
        res['xl'] = xl
        res['xr'] = xr
        res['fl'] = fl
        res['fr'] = fr

        res['nit'] = np.maximum(res['nit'][:n], res['nit'][n:])
        res['nfev'] = res['nfev'][:n] + res['nfev'][n:]
        # If the status on one side is zero, the status is zero. In any case,
        # report the status from one side only.
        res['status'] = np.choose(sa == 0, (sb, sa))
        res['success'] = (res['status'] == 0)

        del res['x']
        del res['f']
        del res['x_last']
        del res['f_last']

        return shape[:-1]

    return eim._loop(work, callback, shape, maxiter, func, args, dtype,
                     pre_func_eval, post_func_eval, check_termination,
                     post_termination_check, customize_result, res_work_pairs)
