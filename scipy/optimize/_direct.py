from __future__ import annotations
from typing import (
    Any, Callable, Iterable, Optional, Tuple, TYPE_CHECKING, Union
)

import numpy as np
from scipy.optimize import OptimizeResult
from ._constraints import old_bound_to_new, Bounds
from ._directmodule import direct as _direct  # type: ignore

if TYPE_CHECKING:
    import numpy.typing as npt
    from _typeshed import NoneType

__all__ = ['direct']

ERROR_MESSAGES = (
    "Number of function evaluations done is larger than maxfun={}",
    "Number of iterations is larger than maxiter={}",
    "u[i] < l[i] for some i",
    "maxfun is too large",
    "Initialization failed",
    "There was an error in the creation of the sample points",
    "An error occured while the function was sampled",
    "Maximum number of levels has been reached.",
    "Forced stop",
    "Invalid arguments",
    "Out of memory",
)

SUCCESS_MESSAGES = (
    ("The best function value found is within a relative error={}"
     "of the (known) global optimum f_min"),
    ("The volume of the hyperrectangle containing the lowest function value "
     "found is below vol_tol={}"),
    ("The side length of the hyperrectangle containing the lowest function "
     "value found is below len_tol={}"),
)


def direct(
    func: Callable[[npt.ArrayLike, Tuple[Any]], float],
    bounds: Union[Iterable, Bounds],
    *,
    args: tuple = (),
    disp: bool = False,
    eps: float = 1e-4,
    maxfun: int = 20000,
    maxiter: int = 6000,
    locally_biased: bool = True,
    f_min: float = -np.inf,
    f_min_rtol: float = 0.0001,
    vol_tol: float = 1e-16,
    len_tol: float = 1e-8,
    callback: Optional[Callable[[npt.ArrayLike], NoneType]] = None
) -> OptimizeResult:
    r"""
    Solve an optimization problem using the DIRECT
    (Dividing Rectangles) algorithm.

    The algorithm is due to Jones et al. [1]_.

    Parameters
    ----------
    func : callable
        The objective function to be minimized.
        ``func(x, *args) -> float``
        where ``x`` is an 1-D array with shape (n,) and ``args`` is a tuple of
        the fixed parameters needed to completely specify the function.
    bounds : sequence or `Bounds`
        Bounds for variables. There are two ways to specify the bounds:
        1. Instance of `Bounds` class.
        2. ``(min, max)`` pairs for each element in ``x``, defining the finite
        lower and upper bounds for the optimizing argument of `func`. It is
        required to have ``len(bounds) == len(x)``. ``len(bounds)`` is used
        to determine the number of parameters in ``x``.
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.
    disp : bool, optional
        If ``True``, print logging information about the optimization process.
    eps : float, optional
        Ensures sufficient decrease in function value when a new potentially
        optimal hyperrectangle is chosen. Default is 0.0001.
    maxfun : int, optional
        Approximate upper bound on objective function evaluations.
        Default is 20000.
    maxiter : int, optional
        Maximum number of iterations. Default is 6000.
    locally_biased : bool, optional
        If `True` (default), use the locally biased variant [2] of the
        algorithm known as DIRECT_L. If `False`, use the original unbiased
        DIRECT algorithm [1]. For hard problems with many local minima,
        `False` is recommended.
    f_min : float, optional
        Function value of the global optimum. Set this value only if the
        global optimum is known. Default is ``-np.inf``, so that this
        termination criterion is deactivated.
    f_min_rtol : float, optional
        Terminate the optimization once the relative error
        `(f - f_min)/f_min` between the current best minimum `f` and
        the supplied global minimum `f_min` is smaller than `f_min_tol`.
        This parameter is only used when `f_min` is also set.
        Default is 0.0001.
    vol_tol : float, optional
        Terminate the optimization once the volume of the hyperrectangle
        containing the lowest function value is smaller than `vol_tol`
        of the complete search space. Must lie between 0 and 1.
        Default is 1e-16.
    len_tol : float, optional
        Terminate the optimization once the maximal side length of the
        hyperrectangle containing the lowest function value is smaller than
        `len_tol`. Must lie between 0 and 1. Default is 1e-8.
    callback : callable, optional
        A callback function with signature ``callback(xk)`` where ``xk``
        represents the best function value found so far.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.

    Notes
    -----

    DIRECT is a deterministic global
    optimization algorithm capable of minimizing a black box function with
    its variables subject to lower and upper bound constrains by sampling
    potential solutions in the search space. The algorithm starts by
    normalising the search space to an n-dimensional unit hypercube.
    It samples the function at the center of this hypercube and at 2n
    (n is the number of variables) more points, 2 in each coordinate
    direction. Using these function values, DIRECT then divides the
    domain into hyperrectangles, each having exactly one of the sampling
    points as its center. In each iteration, DIRECT chooses, using the eps
    parameter which defaults to 1e-4, some of the existing hyperrectangles
    to be further divided. This division process continues until either the
    maximum number of iterations or maximum function evaluations allowed
    are exceeded, or the volume or the side length of the hyperrectangle
    containing the minimal value found so far is lower than a certain
    tolerance. If `f_min` is specified, the optimization will stop once
    this function value is reached within a relative tolerance. The locally
    biased variant of DIRECT (originally called DIRECT_L) [2]_ is used by
    default. It makes the search more locally biased and more efficient for
    cases with only a few local minima.

    This code is based on the DIRECT 2.0.4 Fortran code by Gablonsky et al. at
    https://ctk.math.ncsu.edu/SOFTWARE/DIRECTv204.tar.gz
    The C version was initially converted via f2c and then cleaned up and
    reorganized by Steven G. Johnson, August 2007. This method wraps
    the C implementation.

    .. versionadded:: 1.9.0

    References
    ----------
    .. [1] Jones, D.R., Perttunen, C.D. & Stuckman, B.E. Lipschitzian
           optimization without the Lipschitz constant. J Optim Theory Appl
           79, 157-181 (1993)
    .. [2] Gablonsky, J., Kelley, C. A Locally-Biased form of the DIRECT
           Algorithm. Journal of Global Optimization 21, 27-37 (2001).

    Examples
    ----------
    The following example is a 2-D problem with four local minima: minimizing
    the Styblinski-Tang function
    (https://en.wikipedia.org/wiki/Test_functions_for_optimization).

    >>> from scipy.optimize import direct, Bounds
    >>> def styblinski_tang(pos):
    ...     x, y = pos
    ...     return 0.5 * (x**4 - 16 * x**2 + 5 * x + y**4 - 16 * y**2 + 5 * y)
    >>> bounds = Bounds([-4., -4.], [4., 4.])
    >>> result = direct(styblinski_tang, bounds)
    >>> result.x, result.fun, result.nfev
    array([-2.90362242, -2.90362242]), -78.33233113735979, 20003

    The correct global minimum was found but with a huge number of function
    evaluations (20003). Loosening the termination criteria can be used to stop
    DIRECT earlier.

    >>> from scipy.optimize import direct, Bounds
    >>> def styblinski_tang(pos):
    ...     x, y = pos
    ...     return 0.5 * (x**4 - 16 * x**2 + 5 * x + y**4 - 16 * y**2 + 5 * y)
    >>> bounds = Bounds([-4., -4.], [4., 4.])
    >>> result = direct(styblinski_tang, bounds, vol_tol = 1e-7)
    >>> result.x, result.fun, result.nfev
    array([-2.90321597, -2.9044353 ]), -78.33231560853986, 1113
    """

    if not isinstance(bounds, Bounds):
        lb, ub = old_bound_to_new(bounds)
        bounds = Bounds(lb, ub)

    lb = np.ascontiguousarray(bounds.lb)
    ub = np.ascontiguousarray(bounds.ub)

    def _func_wrap(x, *args):
        x = np.asarray(x)
        f = func(x, *args)
        return f

    x, fun, ret_code, nfev, nit = _direct(
        _func_wrap,
        np.asarray(lb), np.asarray(ub),
        args,
        disp, eps, maxfun, maxiter,
        locally_biased,
        f_min, f_min_rtol,
        vol_tol, len_tol, callback
    )

    format_val = (maxfun, maxiter, f_min_rtol, vol_tol, len_tol)
    if ret_code > 2:
        message = SUCCESS_MESSAGES[ret_code - 3].format(
                    format_val[ret_code - 1])
    elif 0 < ret_code <= 2:
        message = ERROR_MESSAGES[ret_code - 1].format(format_val[ret_code - 1])
    elif 0 > ret_code > -100:
        message = ERROR_MESSAGES[abs(ret_code) + 1]
    else:
        message = ERROR_MESSAGES[ret_code + 99]

    return OptimizeResult(x=np.asarray(x), fun=fun, status=ret_code,
                          success=ret_code > 2, message=message,
                          nfev=nfev, nit=nit)
