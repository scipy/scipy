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
    ("The best function value found is within {} percent "
     "of the (known) global optimum f_min"),
    ("The volume of the hyper-rectangle with best function value found "
     "is below vol_per={}"),
    ("The measure of the hyper-rectangle with best function value found "
     "is below sigma_per={}"),
)


def direct(
    func: Callable[[npt.ArrayLike, Tuple[Any]], float],
    bounds: Union[Iterable, Bounds],
    *,
    args: tuple = (),
    disp: bool = False,
    iatol: float = 1e-4,
    maxfun: int = 20000,
    maxiter: int = 6000,
    locally_biased: bool = True,
    f_min: float = -np.inf,
    f_min_per: float = 0.01,
    vol_per: float = -1.0,
    sigma_per: float = -1.0,
    callback: Optional[Callable[[npt.ArrayLike], NoneType]] = None
) -> OptimizeResult:
    r"""
    Solve an optimization problem using the DIRECT
    (Dividing Rectangles) algorithm.

    The algorithm is due to Jones et al. [1]_.

    Parameters
    ----------
    func: callable
        The objective function to be minimized.
        ``func(x, \*args) -> float``
        where ``x`` is an 1-D array with shape (n,) and args is a tuple of
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
        Prints the evaluated `func` at every iteration.
    iatol : float, optional
        Ensures sufficient decrease in function value when a new potentially
        optimal hyperrectangle is chosen. Default is 0.0001.
    maxfun : int, optional
        Approximate upper bound on objective function evaluations.
        Default is 2000.
    maxiter : int, optional
        Maximum number of iterations. Default is 6000.
    locally_biased : bool, optional
        Whether to use the original [1]_ or modified [2]_ DIRECT algorithm.
        Default is True. Possible values:

        * ``locally_biased=False`` - use the original DIRECT algorithm
        * ``locally_biased=True`` - use the modified DIRECT-l algorithm
    f_min : float, optional
        Function value of the global optimum.
        By default it is a negative value whose absolute value is very large.
        Set this value only if the global optimum is known.
    f_min_per : float, optional
        Terminate the optimization when the percent error satisfies:

        .. math::

            100*(f_{min} - f_{global})/\max(1, |f_{global}|) \leq f_{glper}

        Default is 0.01.
    vol_per : float, optional
        Terminate the optimization once the volume of a hyperrectangle is less
        than vol_per percent of the original hyper-rectangle.
        Default is -1.
    sigma_per : float, optional
        Terminate the optimization once the measure of the hyperrectangle is
        less than this argument. Default is -1.
    callback : callable, `callback(xk)`, optional
        A function to follow the progress of the minimization. ``xk`` is
        the best solution found so far.

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

    DIRECT is a deterministic
    optimization algorithm capable of minimizing a black box function with
    its variables subject to lower and upper bound constrains by sampling
    potential solutions in the search space. The algorithm starts by
    normalising the hyperrectangle (set of all possible values that can
    be taken by the input variables subject to the bound constraints) to
    an n-dimensional unit hypercube. It samples the function at the
    center of this hypercube and at 2n (n is the number of variables)
    more points, 2 in each coordinate direction. Using these function
    values, DIRECT then divides the domain into hyperrectangles, each
    having exactly one of the sampling points as its center. In each
    iteration, DIRECT chooses, using the epsilon parameter, which defaults
    to 1e-4, some of the existing hyperrectangles to be further divided.
    This division process continues until the maximum iterations or
    maximum function evaluations allowed are exceeded, or the function
    value is within the desired percentage error of the global minimum
    (if known). The locally biased variant of DIRECT
    (originally called DIRECT_L) [2]_ is used by default. It makes the
    search more locally biased and more efficient for cases with only a
    few local minima.

    This code is based on the DIRECT 2.0.4 Fortran code by Gablonsky et al. at
    http://www4.ncsu.edu/~ctk/SOFTWARE/DIRECTv204.tar.gz
    The C version was initially converted via f2c and then cleaned up and
    reorganized by Steven G. Johnson, August 2007. And this method wraps
    the C implementation.

    .. versionadded:: 1.8.0

    References
    ----------
    .. [1] Jones, D.R., Perttunen, C.D. & Stuckman, B.E. Lipschitzian
           optimization without the Lipschitz constant. J Optim Theory Appl
           79, 157-181 (1993)
    .. [2] Gablonsky, J., Kelley, C. A Locally-Biased form of the DIRECT
           Algorithm. Journal of Global Optimization 21, 27-37 (2001).
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
        disp, iatol, maxfun, maxiter,
        locally_biased,
        f_min, f_min_per,
        vol_per, sigma_per, callback
    )

    format_val = (maxfun, maxiter, f_min_per, vol_per, vol_per)
    if ret_code > 2:
        message = SUCCESS_MESSAGES[ret_code - 1].format(
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
