import numpy as np
from ._directmodule import direct # type: ignore
from .optimize import OptimizeResult

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


def _minimize_direct(func, bounds, args=(), disp=False,
                     eps=1e-4, maxfun=20000, maxiter=6000,
                     locally_biased=True, f_min=-np.inf, f_min_per=0.01,
                     vol_per=-1.0, sigma_per=-1.0):
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
        fun(x, *args) -> float
        where x is an 1-D array with shape (n,) and args is a tuple of
        the fixed parameters needed to completely specify the function.
    bounds : Bounds
        lower bounds and upper bounds for each element in ``x``, defining
        the bounds on that parameter.
    args : tuple, optional
        Extra arguments passed to the objective function.
    eps : float, optional
        Ensures sufficient decrease in function value when a new potentially
        optimal interval is chosen.
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
        Terminate the optimization once the measure of the
        hyper-rectangle is less than this argument.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.

    References
    ----------
    .. [1] Jones, D.R., Perttunen, C.D. & Stuckman, B.E. Lipschitzian
           optimization without the Lipschitz constant. J Optim Theory Appl
           79, 157-181 (1993)
    .. [2] Jorg Maximilian Xaver Gablonsky and Carl Timothy Kelley. 2001.
           Modifications of the direct algorithm. Ph.D. Dissertation.
    """
    lb = np.ascontiguousarray(bounds.lb)
    ub = np.ascontiguousarray(bounds.ub)

    def _func_wrap(x, *args):
        f = func(x, *args)
        if np.isnan(f):
            return np.nan
        else:
            return f

    x, fun, ret_code, nfev, nit = direct(_func_wrap,
                                         np.asarray(lb), np.asarray(ub),
                                         args,
                                         disp, eps, maxfun, maxiter,
                                         locally_biased,
                                         f_min, f_min_per,
                                         vol_per, sigma_per)

    format_val = (maxfun, maxiter, f_min_per, vol_per, vol_per)
    if ret_code > 2:
        message = SUCCESS_MESSAGES[ret_code - 1].format(
                    format_val[ret_code - 1])
    elif ret_code > 0 and ret_code <= 2:
        message = ERROR_MESSAGES[ret_code - 1].format(format_val[ret_code - 1])
    elif ret_code < 0 and ret_code > -100:
        message = ERROR_MESSAGES[abs(ret_code) + 1]
    else:
        message = ERROR_MESSAGES[ret_code + 99]

    return OptimizeResult(x=np.asarray(x), fun=fun, status=ret_code,
                          success=ret_code > 2, message=message,
                          nfev=nfev, nit=nit)
