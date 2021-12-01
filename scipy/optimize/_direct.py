import numpy as np
from scipy.optimize import OptimizeResult
from ._constraints import old_bound_to_new, Bounds
from ._directmodule import direct as _direct  # type: ignore

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
     "of the (known) global optimum"),
    ("The volume of the hyper-rectangle with best function value found "
     "is below vol_per={}"),
    ("The measure of the hyper-rectangle with best function value found "
     "is below sigma_per={}"),
)


def direct(func, bounds, *, args=(), disp=False,
           iatol=1e-4, maxfun=20000, maxiter=6000,
           locally_biased=True, f_min=-np.inf, f_min_per=0.01,
           vol_per=-1.0, sigma_per=-1.0, callback=None):
    r"""
    Solve an optimization problem using the DIRECT
    (Dividing Rectangles) algorithm.

    It can be used to solve general nonlinear programming problems of the form:

    .. math::

           \min_ {x \in R^n} f(x)

    subject to

    .. math::

           x_L \leq  x  \leq x_U

    Where :math:`x` are the optimization variables (with upper and lower
    bounds), :math:`f(x)` is the objective function.

    Parameters
    ----------
    fun: callable
        The objective function to be minimized.
        ``fun(x, \*args) -> float``
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
        Extra arguments passed to the objective function.
    iatol : float, optional
        Ensures sufficient decrease in function value when a new potentially
        optimal hyper-rectangle is chosen.
    maxfun : int, optional
        Approximate upper bound on objective function evaluations.
    maxiter : int, optional
        Maximum number of iterations.
    locally_biased : bool, optional
        Whether to use the original [1]_ or modified [2]_ DIRECT algorithm.
        Possible values:

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
    vol_per : float, optional
        Terminate the optimization once the volume of a hyper-rectangle is less
        than vol_per percent of the original hyper-rectangle.
    sigma_per : float, optional
        Terminate the optimization once the measure of the hyper-rectangle is
        less than this argument.

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
    optimization algorithm capable of minimizing black box function with
    its variables subject to lower and upper bound constrains by sampling
    potential solutions in the search space. The algorithm starts with
    mapping the hyperrectangle (set of all possible values that can be taken
    by the input variables subject to the bound constraints) n-dimensional
    unit hypercube. It samples the function at the center of this hypercube
    and at 2n (n is the number of variables) more points, 2 in each coordinate
    direction. Using these function values, DIRECT then divides the domain
    into hyperrectangles, each having exactly one of the sampling points as
    its center. In each iteration, DIRECT chooses, using the epsilon parameter,
    by default 1e-4, some of the existing hyperrectangles to be further divided
    . This division process continues until the maximum iterations or maximum
    function evaluations allowed are exceeded, or the function value is within
    the desired percentage error of the global minimum (if known). The improved
    version of DIRECT algorithm [1]_ is biased towards local search making it
    effective for functions without too many local minima. This method wraps
    the C implementation of the original and improved algorithms.

    .. versionadded:: 1.8.0

    References
    ----------
    .. [1] Jones, D.R., Perttunen, C.D. & Stuckman, B.E. Lipschitzian
           optimization without the Lipschitz constant. J Optim Theory Appl
           79, 157-181 (1993)
    .. [2] Jorg Maximilian Xaver Gablonsky and Carl Timothy Kelley. 2001.
           Modifications of the direct algorithm. Ph.D. Dissertation.
    """
    if not isinstance(bounds, Bounds):
        lb, ub = old_bound_to_new(bounds)
        bounds = Bounds(lb, ub)

    lb = np.ascontiguousarray(bounds.lb)
    ub = np.ascontiguousarray(bounds.ub)

    def _func_wrap(x, *args):
        x = np.asarray(x)
        f = func(x, *args)
        if np.isnan(f):
            return np.nan
        else:
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
