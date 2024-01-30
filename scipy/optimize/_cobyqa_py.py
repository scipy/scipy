import numpy as np

from ._optimize import _check_unknown_options

__all__ = []


def _minimize_cobyqa(fun, x0, args=(), bounds=None, constraints=(),
                     callback=None, disp=False, maxfev=1000, maxiter=1000,
                     target=-np.inf, feasibility_tol=1e-8, radius_init=1.0,
                     tol=1e-6, scale=False, **unknown_options):
    """
    Minimize a scalar function of one or more variables using the
    Constrained Optimization BY Quadratic Approximations (COBYQA) algorithm.

    Options
    -------
    disp : bool
        Set to True to print information about the optimization procedure.
    maxfev : int
        Maximum number of function evaluations.
    maxiter : int
        Maximum number of iterations.
    target : float
        Target value for the objective function. The optimization procedure is
        terminated when the objective function value of a nearly feasible point
        is less than or equal to this target.
    feasibility_tol : float
        Tolerance for the constraint violation.
    radius_init : float
        Initial trust-region radius. Typically, this value should be in the
        order of one tenth of the greatest expected change to the variables.
    tol : float
        Final trust-region radius. It should indicate the accuracy required in
        the final values of the variables.
    scale : bool
        Whether to scale the variables according to the bounds.
    """
    from .._lib.cobyqa import minimize  # import here to avoid circular imports

    _check_unknown_options(unknown_options)
    options = {
        'disp': bool(disp),
        'maxfev': int(maxfev),
        'maxiter': int(maxiter),
        'target': float(target),
        'feasibility_tol': float(feasibility_tol),
        'radius_init': float(radius_init),
        'radius_final': float(tol),
        'scale': bool(scale),
    }
    return minimize(fun, x0, args, bounds, constraints, callback, options)
