from typing import (  # noqa: UP035
    Any, Callable, Iterable
)

import numpy as np
from scipy.optimize import OptimizeResult
from ._constraints import old_bound_to_new, Bounds
from ._biteopt import minimize as _minimize  # type: ignore
from scipy._lib._util import _validate_int

__all__ = ['biteopt']


def biteopt(
    func: Callable[[np.ndarray, Any], float],
    bounds: Iterable | Bounds,
    *,
    args: tuple = (),
    maxfun: int | None = None,
    depth: int = 1,
    f_min: float | None = None,
    rng: int | np.random.Generator | None = None,
):
    """Find the global minimum of a function using the BiteOpt algorithm.

    Parameters
    ----------
    func : callable
        The objective function to be minimized, ``func(x, *args) -> float``,
        where ``x`` is a 1-D array with shape ``(n,)`` and ``args`` is a tuple
        of fixed parameters.
    bounds : sequence or `Bounds`
        Bounds for variables, specified either as an instance of `Bounds` or as
        ``(min, max)`` pairs for each element in ``x``. Bounds must be finite.
    args : tuple, optional
        Additional fixed parameters passed to the objective function.
    maxfun : int, optional
        Maximum number of objective function evaluations. Default is
        ``1000 * n``, where ``n`` is the number of variables inferred from
        ``bounds``.
    depth : int, optional
        Algorithm depth ``M``. ``1`` selects the plain BiteOpt algorithm, while
        values ``> 1`` select the deeper BiteOptDeep variant. Expected range is
        ``[1, 36]``. Default is 1.
    f_min : float, optional
        Target objective value. The optimization stops early once the best
        objective value found is less than or equal to `f_min`. By default
        (`None`) this criterion is disabled and the full iteration budget is
        used.
    rng : {None, int, `numpy.random.Generator`}, optional
        Controls reproducibility. Passed to `numpy.random.default_rng` to
        derive the integer seed used by BiteOpt's internal PRNG.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as an `OptimizeResult` object.
        Important attributes are: ``x`` the solution array, ``fun`` the value
        of the objective at the solution, ``nfev`` the number of objective
        evaluations performed, ``success`` a boolean flag indicating whether
        the optimizer terminated successfully, and ``message`` describing the
        cause of termination. When `f_min` is given, ``success`` is ``True``
        only if that target was reached. Absent `f_min`, BiteOpt has no
        reliable convergence test and runs its full iteration budget, so a
        completed run reports ``success`` as ``True``.

    Notes
    -----
    BiteOpt is a stochastic, derivative-free global optimizer that maintains a
    small population of candidate solutions and evolves them with an ensemble
    of stochastic move generators governed by an adaptive selection scheme. It
    targets low- to medium-dimensional continuous problems with finite box
    bounds and requires no gradient information. Because the search is
    stochastic, results depend on the random stream; pass `rng` for
    reproducible runs.

    .. versionadded:: 1.19.0

    References
    ----------
    .. [1] Aleksey Vaneev. "BiteOpt - Derivative-Free Global Optimization
           Method (C++)". https://github.com/avaneev/biteopt

    Examples
    --------
    The following example is a 2-D problem with four local minima: minimizing
    the Styblinski-Tang function
    (https://en.wikipedia.org/wiki/Test_functions_for_optimization).

    >>> import numpy as np
    >>> from scipy.optimize import biteopt, Bounds
    >>> def styblinski_tang(pos):
    ...     x, y = pos
    ...     return 0.5 * (x**4 - 16*x**2 + 5*x + y**4 - 16*y**2 + 5*y)
    >>> bounds = Bounds([-4., -4.], [4., 4.])
    >>> result = biteopt(styblinski_tang, bounds)
    >>> result.x, result.fun, result.nfev
    array([-2.90353406, -2.90353401]), -78.3323279095383, 2000  # may vary

    For reproducible results, pass a seed to `rng`:
    >>> rng = np.random.default_rng(1234)
    >>> result = biteopt(styblinski_tang, bounds, rng=rng)
    >>> result.x, result.fun, result.nfev
    array([-2.90353404, -2.90353406]), -78.33233140754282, 2000. # may vary

    To stop the optimization early once a target objective value is reached,
    pass `f_min`:
    >>> result = biteopt(styblinski_tang, bounds, f_min=-70, rng=rng)
    >>> result.x, result.fun, result.nfev
    array([-3.06877089, -3.06877089]), -77.33495993486557, 31  # may vary
    """
    if not callable(func):
        raise TypeError("func must be callable")

    if not isinstance(bounds, Bounds):
        if isinstance(bounds, (list, tuple)):
            if len(bounds) == 0:
                raise ValueError(
                    "bounds must contain at least one finite (min, max) pair"
                )
            lb, ub = old_bound_to_new(bounds)
            bounds = Bounds(lb, ub)
        else:
            raise ValueError(
                "bounds must be a sequence or instance of Bounds class"
            )

    lb = np.ascontiguousarray(bounds.lb, dtype=np.float64)
    ub = np.ascontiguousarray(bounds.ub, dtype=np.float64)

    if lb.shape != ub.shape or lb.ndim != 1 or lb.shape[0] == 0:
        raise ValueError("bounds must contain at least one finite (min, max) pair")
    if not np.all(lb < ub):
        raise ValueError("Bounds are not consistent min < max")
    if np.any(np.isinf(lb)) or np.any(np.isinf(ub)):
        raise ValueError("Bounds must not be inf.")

    depth = _validate_int(depth, "depth", minimum=1)
    if depth > 36:
        raise ValueError("depth must be an integer in [1, 36].")

    if maxfun is not None:
        maxfun = _validate_int(maxfun, "maxfun", minimum=1)

    if f_min is not None:
        f_min = float(f_min)

    # BiteOpt's internal PRNG is driven by the NumPy Generator's bit
    # generator. Keep a reference to ``generator`` alive for the duration of
    # the call so its capsule pointer stays valid.
    generator = np.random.default_rng(rng)

    def func_wrap(x):
        fx = func(x, *args)
        if not np.isscalar(fx):
            _dt = getattr(fx, "dtype", np.dtype(np.float64))
            try:
                fx = _dt.type(np.asarray(fx).item())
            except (TypeError, ValueError) as e:
                raise ValueError(
                    "The user-provided objective function "
                    "must return a scalar value."
                ) from e
        return fx

    if maxfun is None:
        maxfun = 1000 * len(lb)

    # naming convention: SciPy's "maxfun" is BiteOpt's "iter"
    # BiteOpt internally scales iter by sqrt(depth): useiter = int(iter * sqrt(depth)).
    # Pre-divide so the actual evaluation count never exceeds maxfun.
    # floor(maxfun / sqrt(depth)) guarantees no overshoot; then probe floor+1
    # to squeeze in one more iteration if it still fits within the budget.
    _sqrt = np.sqrt(depth)
    _iter = int(maxfun / _sqrt)
    if int((_iter + 1) * _sqrt) <= maxfun:
        _iter += 1
    _iter = max(1, _iter)

    result = _minimize(
        func_wrap,
        lb,
        ub,
        _iter,
        int(depth),
        1, # only one attempt, the best result is returned
        generator.bit_generator.capsule,
        f_min,
    )

    # BiteOpt is a fixed-budget stochastic optimizer. Its only intrinsic success
    # criterion is the early-stop target `f_min`, which it reaches exactly when
    # the best value satisfies ``fun <= f_min``. Absent `f_min` there is no
    # reliable convergence test (biteopt typically keeps making small improving
    # steps until the iteration budget is exhausted), so a completed run is
    # reported as successful.
    if f_min is None:
        success = True
        message = "Maximum number of iterations reached."
    elif result["fun"] <= f_min:
        success = True
        message = "Optimization terminated successfully: f_min reached."
    else:
        success = False
        message = "Maximum number of iterations reached; f_min not reached."

    return OptimizeResult(
        x=np.asarray(result["x"]),
        fun=result["fun"],
        nfev=result["nfev"],
        success=success,
        message=message,
    )
