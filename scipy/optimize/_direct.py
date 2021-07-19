r"""
A python wrapper to the DIRECT algorithm implemented in Fortran
===============================================================

DIRECT is a method to solve global bound constraint optimization problems and
was originally developed by D. R. Jones, C. D. Perttunen and B. E. Stuckmann.
It is designed to find **global** solutions of mathematical optimization 
problems of the from,

.. math::

       \min_ {x \in R^n} f(x)

subject to

.. math::

       x_L \leq  x  \leq x_U

Where :math:`x` are the optimization variables (with upper and lower
bounds), :math:`f(x)` is the objective function.

The DIRECT package uses the Fortran implementation of DIRECT written by
Joerg.M.Gablonsky, DIRECT Version 2.0.4. More information on the DIRECT
algorithm can be found in Gablonsky's `thesis 
<http://repository.lib.ncsu.edu/ir/bitstream/1840.16/3920/1/etd.pdf>`_.

Authors:

Andreas Mayer <andimscience@gmail.com>, Amit Aides <amitibo@tx.technion.ac.il>
"""

import numpy as np
from ._directmodule import direct
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
     "is below volper={}"),
    ("The measure of the hyper-rectangle with best function value found "
     "is below sigmaper={}"),
)

def _minimize_direct(func, bounds, *args, disp=False,
                     eps=1e-4, maxfun=2000, maxiter=6000, 
                     locally_biased=False, fglobal=-1e100, fglper=0.01, 
                     volper=-1.0, sigmaper=-1.0):
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
        where x is an 1-D array with shape (n,) and args is a tuple of the fixed parameters needed to completely specify the function.
    bounds : Bounds
        lower bounds and upper bounds for each element in ``x``, defining
        the bounds on that parameter.
    eps : float, optional
        Ensures sufficient decrease in function value when a new potentially
        optimal interval is chosen.
    maxfun : int, optional
        Approximate upper bound on objective function evaluations.
    maxiter : int, optional
        Maximum number of iterations.
    locally_biased : bool, optional
        Whether to use the original or modified DIRECT algorithm.
        Possible values:

        * ``locally_biased=False`` - use the original DIRECT algorithm
        * ``locally_biased=True`` - use the modified DIRECT-l algorithm
    fglobal : float, optional
        Function value of the global optimum.
        By default it is a negative value whose absolute value is very large.
        Set this value only if the global optimum is known.
    fglper : float, optional
        Terminate the optimization when the percent error satisfies:

        .. math::

            100*(f_{min} - f_{global})/\max(1, |f_{global}|) \leq f_{glper}
    volper : float, optional
        Terminate the optimization once the volume of a hyperrectangle is less
        than volper percent of the original hyperrectangel.
    sigma_per : float, optional
        Terminate the optimization once the measure of the 
        hyperrectangle is less than this argument.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.
    """
    l = np.ascontiguousarray(bounds.lb)
    u = np.ascontiguousarray(bounds.ub)

    def func_wrap(x, *args):
        try:
            return func(x, *args)
        except Exception:
            return np.nan

    #
    # Call the DIRECT algorithm
    #
    x, fun, ret_code, nfev, nit = direct(func_wrap, np.asarray(l), np.asarray(u), args,
        disp, eps, maxfun, maxiter,
        locally_biased,
        fglobal, fglper,
        volper, sigmaper)

    format_val = (maxfun, maxiter, fglper, volper, volper)
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
            success=ret_code > 2, message=message, nfev=nfev, nit=nit)
