from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc
from .constraints import (NonlinearConstraint,
                          LinearConstraint,
                          BoxConstraint,
                          CanonicalConstraint,
                          concatenate_canonical_constraints,
                          empty_canonical_constraint,
                          generate_lagrangian_hessian)
from ._tr_interior_point import tr_interior_point
from ._equality_constrained_sqp import equality_constrained_sqp

TERMINATION_MESSAGES = {
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`xtol` termination condition is satisfied.",
    3: "`callback` function requested termination"
}


def minimize_constrained(fun, x0, grad, hess=None, constraints=(),
                         method=None, xtol=1e-8, gtol=1e-8,
                         options={}, callback=None, max_iter=1000,
                         sparse_jacobian=None):
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
    constraints : Constraint or List of Constraint's
        A single constraint or a list of constraints.
        Available constraints are:

            - BoxConstraint
            - LinearConstraint
            - NonlinearConstraint

    xtol : float, optional
        Tolerance for termination by the change of the independent variable.
        Only for 'equality-constrained-sqp'.
    gtol : float, optional
        Tolerance for termination by the norm of the lagrangian gradient.
    options : dict, optional
        A dictionary of solver options.

            - `return_all`

    callback : callable, optional
        Called after each iteration:

            callback(OptimizeResult state) -> bool

        and if callback returns true the algorithm execution is terminated.
    max_iter : int, optional
        Maximum number of algorithm iterations.
    sparse_jacobian : {boolean, None}
        The algorithm uses a sparse representation of the Jacobian if True
        and a dense representation if False. When ``None`` the algorithm
        uses the more convenient option (which is not always the more effective
        one).

    Returns
    -------
    result : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.
    """
    # Put ``constraints`` in list format
    if isinstance(constraints, (NonlinearConstraint,
                                LinearConstraint,
                                BoxConstraint,
                                CanonicalConstraint)):
        constraints = [constraints]
    if isinstance(constraints, (list, tuple, np.array)):
        # Converts all constraints to canonical format
        constraints_list = []
        for constr in constraints:
            if not isinstance(constr, (NonlinearConstraint,
                                       LinearConstraint,
                                       BoxConstraint,
                                       CanonicalConstraint)):
                raise ValueError("Unknown Constraint type")
            elif isinstance(constr, CanonicalConstraint):
                constraints_list += [constr]
            else:
                constraints_list += [constr.to_canonical(sparse_jacobian)]
        # Concatenate constraints
        if len(constraints_list) == 0:
            constr = empty_canonical_constraint()
        elif len(constraints_list) == 1:
            constr = constraints_list[0]
        else:
            constr = concatenate_canonical_constraints(constraints_list,
                                                       sparse_jacobian)
    else:
        raise ValueError("Unknown Constraint type")

    # Generate lagrangian hess function
    lagr_hess = generate_lagrangian_hessian(constr, hess)

    # Choose appropriate method
    if method is None:
        if constr.n_ineq == 0:
            method = 'equality_constrained_sqp'
        else:
            method = 'tr_interior_point'

    # Define stop criteria
    if method == 'equality_constrained_sqp':
        def stop_criteria(state):
            state.status = None
            if (callback is not None) and callback(state):
                state.status = 3
            elif state.optimality < gtol and state.constr_violation < gtol:
                state.status = 1
            elif state.trust_radius < xtol:
                state.status = 2
            elif state.niter > max_iter:
                state.status = 0
            return state.status in (0, 1, 2, 3)
    elif method == 'tr_interior_point':
        def stop_criteria(state):
            state.status = None
            if (callback is not None) and callback(state):
                state.status = 3
            elif state.optimality < gtol and state.constr_violation < gtol:
                state.status = 1
            elif state.niter > max_iter:
                state.status = 0
            return state.status in (0, 1, 2, 3)

    # Call inferior function to do the optimization
    if method == 'equality_constrained_sqp':
        result = equality_constrained_sqp(
            fun, grad, lagr_hess,
            constr.constr_eq,  constr.jac_eq,
            x0, stop_criteria, **options)
    elif method == 'tr_interior_point':
        result = tr_interior_point(
            fun, grad, lagr_hess,
            constr.n_ineq, constr.constr_ineq,
            constr.jac_ineq, constr.n_eq,
            constr.constr_eq, constr.jac_eq,
            x0, stop_criteria, None, sparse_jacobian,
            **options)

    result.message = TERMINATION_MESSAGES[result.status]
    result.success = result.status > 0

    return result
