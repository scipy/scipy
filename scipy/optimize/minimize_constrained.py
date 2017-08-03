from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc
from .constraints import (NonlinearConstraint,
                          LinearConstraint,
                          BoxConstraint,
                          concatenate_canonical_constraints)


def minimize_constrained(fun, x0, grad, hess=None, method=None, 
                         constraints=()):
    """Minimize scalar function subject to constraints.

    Parameters
    ----------
    fun : callable
        The objective function to be minimized.

            fun(x) -> float

        where x is an array with shape (n,).
    grad : callable
        Gradient of the objective function:

            grad(x) -> array_like, shape (n,)

        where x is an array with shape (n,).
    hess : {callable, None}, optional
        Hessian of the objective function:

            hess(x) -> LinearOperator (or sparse matrix or ndarray), shape (n, n)

        where x is an array with shape (n,).
    x0 : ndarray
        Initial guess. ``len(x0)`` is the dimensionality of the minimization
        problem.
    method : {str, None}, optional
        Type of solver. Should be on of:

            - 'equality-constrained-sqp'
            - 'tr-interior-point'

        When ``None`` the more appropriate method is choosen.
    constraints : List of constraints
        A list of constraints. Available constraints are:

            - BoxConstraint
            - LinearConstraint
            - NonlinearConstraint

    xtol : float, optional
        Tolerance for termination by the change of the independent variable.
    gtol : float, optional
        Tolerance for termination by the norm of the lagrangian gradient.
    callback : callable, optional
        Called after each iteration:

            callback(OptimizeResult state) -> bool

        and if callback returns true the algorithm execution is interupted.

    Returns
    -------
    result : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.
    """
    pass
