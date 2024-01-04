import numpy as np
from scipy.optimize._optimize import OptimizeResult, _call_callback_maybe_halt

_ESIGNERR = -1
_ECONVERR = -2
_EVALUEERR = -3
_ECALLBACK = -4
_ECONVERGED = 0
_EINPROGRESS = 1


def _elementwise_algorithm_initialize(func, xs, args, complex_ok=False):
    """Initialize abscissa, function, and args arrays for elementwise function

    Parameters
    ----------
    func : callable
        An elementwise function with signature

            func(x: ndarray, *args) -> ndarray

        where each element of ``x`` is a finite real and ``args`` is a tuple,
        which may contain an arbitrary number of arrays that are broadcastable
        with ``x``.
    xs : tuple of arrays
        Finite real abscissa arrays. Must be broadcastable.
    args : tuple, optional
        Additional positional arguments to be passed to `func`.

    Returns
    -------
    xs, fs, args : tuple of arrays
        Broadcasted, writeable, 1D abscissa and function value arrays (or
        NumPy floats, if appropriate). The dtypes of the `xs` and `fs` are
        `xfat`; the dtype of the `args` are unchanged.
    shape : tuple of ints
        Original shape of broadcasted arrays.
    xfat : NumPy dtype
        Result dtype of abscissae, function values, and args determined using
        `np.result_type`, except integer types are promoted to `np.float64`.

    Raises
    ------
    ValueError
        If the result dtype is not that of a real scalar

    Notes
    -----
    Useful for initializing the input of SciPy functions that accept
    an elementwise callable, abscissae, and arguments; e.g.
    `scipy.optimize._chandrupatla`.
    """
    nx = len(xs)

    # Try to preserve `dtype`, but we need to ensure that the arguments are at
    # least floats before passing them into the function; integers can overflow
    # and cause failure.
    # There might be benefit to combining the `xs` into a single array and
    # calling `func` once on the combined array. For now, keep them separate.
    xas = np.broadcast_arrays(*xs, *args)  # broadcast and rename
    xat = np.result_type(*[xa.dtype for xa in xas])
    xat = np.float64 if np.issubdtype(xat, np.integer) else xat
    xs, args = xas[:nx], xas[nx:]
    xs = [x.astype(xat, copy=False)[()] for x in xs]
    fs = [np.asarray(func(x, *args)) for x in xs]
    shape = xs[0].shape

    message = ("The shape of the array returned by `func` must be the same as "
               "the broadcasted shape of `x` and all other `args`.")
    shapes_equal = [f.shape == shape for f in fs]
    if not np.all(shapes_equal):
        raise ValueError(message)

    # These algorithms tend to mix the dtypes of the abscissae and function
    # values, so figure out what the result will be and convert them all to
    # that type from the outset.
    xfat = np.result_type(*([f.dtype for f in fs] + [xat]))
    if not complex_ok and not np.issubdtype(xfat, np.floating):
        raise ValueError("Abscissae and function output must be real numbers.")
    xs = [x.astype(xfat, copy=True)[()] for x in xs]
    fs = [f.astype(xfat, copy=True)[()] for f in fs]

    # To ensure that we can do indexing, we'll work with at least 1d arrays,
    # but remember the appropriate shape of the output.
    xs = [x.ravel() for x in xs]
    fs = [f.ravel() for f in fs]
    args = [arg.flatten() for arg in args]
    return xs, fs, args, shape, xfat


def _elementwise_algorithm_loop(work, callback, shape, maxiter,
                                func, args, dtype, pre_func_eval, post_func_eval,
                                check_termination, post_termination_check,
                                customize_result, res_work_pairs):
    """Main loop of a vectorized scalar optimization algorithm

    Parameters
    ----------
    work : OptimizeResult
        All variables that need to be retained between iterations. Must
        contain attributes `nit`, `nfev`, and `success`
    callback : callable
        User-specified callback function
    shape : tuple of ints
        The shape of all output arrays
    maxiter :
        Maximum number of iterations of the algorithm
    func : callable
        The user-specified callable that is being optimized or solved
    args : tuple
        Additional positional arguments to be passed to `func`.
    dtype : NumPy dtype
        The common dtype of all abscissae and function values
    pre_func_eval : callable
        A function that accepts `work` and returns `x`, the active elements
        of `x` at which `func` will be evaluated. May modify attributes
        of `work` with any algorithmic steps that need to happen
         at the beginning of an iteration, before `func` is evaluated,
    post_func_eval : callable
        A function that accepts `x`, `func(x)`, and `work`. May modify
        attributes of `work` with any algorithmic steps that need to happen
         in the middle of an iteration, after `func` is evaluated but before
         the termination check.
    check_termination : callable
        A function that accepts `work` and returns `stop`, a boolean array
        indicating which of the active elements have met a termination
        condition.
    post_termination_check : callable
        A function that accepts `work`. May modify `work` with any algorithmic
        steps that need to happen after the termination check and before the
        end of the iteration.
    customize_result : callable
        A function that accepts `res` and `shape` and returns `shape`. May
        modify `res` (in-place) according to preferences (e.g. rearrange
        elements between attributes) and modify `shape` if needed.
    res_work_pairs : list of (str, str)
        Identifies correspondence between attributes of `res` and attributes
        of `work`; i.e., attributes of active elements of `work` will be
        copied to the appropriate indices of `res` when appropriate. The order
        determines the order in which OptimizeResult attributes will be
        pretty-printed.

    Returns
    -------
    res : OptimizeResult
        The final result object

    Notes
    -----
    Besides providing structure, this framework provides several important
    services for a vectorized optimization algorithm.

    - It handles common tasks involving iteration count, function evaluation
      count, a user-specified callback, and associated termination conditions.
    - It compresses the attributes of `work` to eliminate unnecessary
      computation on elements that have already converged.

    """
    cb_terminate = False

    # Initialize the result object and active element index array
    n_elements = int(np.prod(shape))
    active = np.arange(n_elements)  # in-progress element indices
    res_dict = {i: np.zeros(n_elements, dtype=dtype) for i, j in res_work_pairs}
    res_dict['success'] = np.zeros(n_elements, dtype=bool)
    res_dict['status'] = np.full(n_elements, _EINPROGRESS)
    res_dict['nit'] = np.zeros(n_elements, dtype=int)
    res_dict['nfev'] = np.zeros(n_elements, dtype=int)
    res = OptimizeResult(res_dict)
    work.args = args

    active = _elementwise_algorithm_check_termination(
        work, res, res_work_pairs, active, check_termination)

    if callback is not None:
        temp = _elementwise_algorithm_prepare_result(
            work, res, res_work_pairs, active, shape, customize_result)
        if _call_callback_maybe_halt(callback, temp):
            cb_terminate = True

    while work.nit < maxiter and active.size and not cb_terminate and n_elements:
        x = pre_func_eval(work)

        if work.args and work.args[0].ndim != x.ndim:
            # `x` always starts as 1D. If the SciPy function that uses
            # _elementwise_algorithm_loop added dimensions to `x`, we need to
            # add them to the elements of `args`.
            dims = np.arange(x.ndim, dtype=np.int64)
            work.args = [np.expand_dims(arg, tuple(dims[arg.ndim:]))
                         for arg in work.args]

        f = func(x, *work.args)
        f = np.asarray(f, dtype=dtype)
        work.nfev += 1 if x.ndim == 1 else x.shape[-1]

        post_func_eval(x, f, work)

        work.nit += 1
        active = _elementwise_algorithm_check_termination(
            work, res, res_work_pairs, active, check_termination)

        if callback is not None:
            temp = _elementwise_algorithm_prepare_result(
                work, res, res_work_pairs, active, shape, customize_result)
            if _call_callback_maybe_halt(callback, temp):
                cb_terminate = True
                break
        if active.size == 0:
            break

        post_termination_check(work)

    work.status[:] = _ECALLBACK if cb_terminate else _ECONVERR
    return _elementwise_algorithm_prepare_result(
        work, res, res_work_pairs, active, shape, customize_result)


def _elementwise_algorithm_check_termination(work, res, res_work_pairs, active,
                                             check_termination):
    # Checks termination conditions, updates elements of `res` with
    # corresponding elements of `work`, and compresses `work`.

    stop = check_termination(work)

    if np.any(stop):
        # update the active elements of the result object with the active
        # elements for which a termination condition has been met
        _elementwise_algorithm_update_active(work, res, res_work_pairs, active,
                                             stop)

        # compress the arrays to avoid unnecessary computation
        proceed = ~stop
        active = active[proceed]
        for key, val in work.items():
            work[key] = val[proceed] if isinstance(val, np.ndarray) else val
        work.args = [arg[proceed] for arg in work.args]

    return active


def _elementwise_algorithm_update_active(work, res, res_work_pairs, active,
                                         mask=None):
    # Update `active` indices of the arrays in result object `res` with the
    # contents of the scalars and arrays in `update_dict`. When provided,
    # `mask` is a boolean array applied both to the arrays in `update_dict`
    # that are to be used and to the arrays in `res` that are to be updated.
    update_dict = {key1: work[key2] for key1, key2 in res_work_pairs}
    update_dict['success'] = work.status == 0

    if mask is not None:
        active_mask = active[mask]
        for key, val in update_dict.items():
            res[key][active_mask] = val[mask] if np.size(val) > 1 else val
    else:
        for key, val in update_dict.items():
            res[key][active] = val


def _elementwise_algorithm_prepare_result(work, res, res_work_pairs, active,
                                          shape, customize_result):
    # Prepare the result object `res` by creating a copy, copying the latest
    # data from work, running the provided result customization function,
    # and reshaping the data to the original shapes.
    res = res.copy()
    _elementwise_algorithm_update_active(work, res, res_work_pairs, active)

    shape = customize_result(res, shape)

    for key, val in res.items():
        res[key] = np.reshape(val, shape)[()]
    res['_order_keys'] = ['success'] + [i for i, j in res_work_pairs]
    return OptimizeResult(**res)
