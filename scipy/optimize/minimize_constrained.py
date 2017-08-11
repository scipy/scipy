from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc
from .constraints import (NonlinearConstraint,
                          LinearConstraint,
                          BoxConstraint,
                          CanonicalConstraint,
                          concatenate_canonical_constraints,
                          empty_canonical_constraint,
                          generate_lagrangian_hessian,
                          reinforce_box_constraints)
from ._tr_interior_point import tr_interior_point
from ._equality_constrained_sqp import equality_constrained_sqp
from warnings import warn
import time
from .optimize import OptimizeResult


TERMINATION_MESSAGES = {
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`xtol` termination condition is satisfied.",
    3: "`callback` function requested termination"
}

class sqp_printer:
    @staticmethod
    def print_header():
        print("|{0:^7}|{1:^7}|{2:^7}|{3:^10}|{4:^10}|{5:^10}|{6:^10}|"
              .format("niter", "f evals", "CG iter", "tr radius", "penalty",
                      "opt", "c viol"))
        s = "-"*6 + ":"
        s2 = ":" + "-"*8 + ":"
        print("|{0:^7}|{1:^7}|{2:^7}|{3:^10}|{4:^10}|{5:^10}|{6:^10}|"
              .format(s, s, s, s2, s2, s2, s2))

    @staticmethod
    def print_problem_iter(niter, nfev, cg_niter, tr_radius,
                           penalty, opt, c_viol):
        print("|{0:>7}|{1:>7}|{2:>7}| {3:^1.2e} | {4:^1.2e} | {5:^1.2e} | {6:^1.2e} |"
              .format(niter, nfev, cg_niter, tr_radius, penalty,
                      opt, c_viol))

    @staticmethod
    def print_footer():
        print("")
        print((7*3 + 10*4 + 8)*"-")
        print("")

class ip_printer:
    @staticmethod
    def print_header():
        print("|{0:^7}|{1:^7}|{2:^7}|{3:^13}|{4:^10}|{5:^10}|{6:^10}|{7:^10}|"
              .format("niter", "f evals", "CG iter", "barrier param",
                      "tr radius", "penalty", "opt", "c viol"))
        s = "-"*6 + ":"
        s2 = ":" + "-"*11 + ":"
        s3 = ":" + "-"*8 + ":"
        print("|{0:^7}|{1:^7}|{2:^7}|{3:^13}|{4:^10}|{5:^10}|{6:^10}|{7:^10}|"
              .format(s, s, s, s2, s3, s3, s3, s3))

    @staticmethod
    def print_problem_iter(niter, nfev, cg_niter, barrier_parameter, tr_radius,
                           penalty, opt, c_viol):
        print("|{0:>7}|{1:>7}|{2:>7}|   {3:^1.2e}  | "
              "{4:^1.2e} | {5:^1.2e} | {6:^1.2e} | {7:^1.2e} |"
              .format(niter, nfev, cg_niter, barrier_parameter,
                      tr_radius, penalty, opt, c_viol))

    @staticmethod
    def print_footer():
        print("")
        print((7*3 + 13 + 10*4 + 9)*"-")
        print("")


def minimize_constrained(fun, x0, grad, hess=None, constraints=(),
                         method=None, xtol=1e-8, gtol=1e-8,
                         options={}, callback=None, max_iter=1000,
                         sparse_jacobian=None, verbose=0):
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
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

            * 0 (default) : work silently.
            * 1 : display a termination report.
            * 2 : display progress during iterations.

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
        for c in constraints:
            if not isinstance(c, (NonlinearConstraint,
                                  LinearConstraint,
                                  BoxConstraint,
                                  CanonicalConstraint)):
                raise ValueError("Unknown Constraint type")
            elif isinstance(c, CanonicalConstraint):
                constraints_list += [c]
            else:
                constraints_list += [c.to_canonical(sparse_jacobian)]
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

    # Compute initial values
    fun0 = fun(x0)
    grad0 = grad(x0)
    constr_eq0 = constr.constr_eq(x0)
    constr_ineq0 = constr.constr_ineq(x0)
    jac_eq0 = constr.jac_eq(x0)
    jac_ineq0 = constr.jac_ineq(x0)

    # Construct OptimizeResult
    state = OptimizeResult(niter=0, nfev=1, ngev=1,
                           ncev=1, njev=1, nhev=0,
                           cg_niter=0, cg_info={})
    # Store values
    return_all = options.get("return_all", False)
    if return_all:
        state.allvecs = []
        state.allmult = []

    # Check initial point
    x0 = np.asarray(x0)
    for c in constraints:
        if np.any(c.feasible_constr):
            if isinstance(c, BoxConstraint):
                x0_new =  reinforce_box_constraints(c, x0)
                if not np.array_equal(x0_new, x0):
                    warn('The initial point was changed in order '
                         +'to stay inside box constraints.')
                    x0 = x0_new

    # Choose appropriate method
    if method is None:
        if constr.n_ineq == 0:
            method = 'equality_constrained_sqp'
        else:
            method = 'tr_interior_point'

    # Define stop criteria
    if method == 'equality_constrained_sqp':
        def stop_criteria(state):
            if verbose >= 2:
                sqp_printer.print_problem_iter(state.niter,
                                               state.nfev,
                                               state.cg_niter,
                                               state.trust_radius,
                                               state.penalty,
                                               state.optimality,
                                               state.constr_violation)
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
            barrier_tol = options.get("barrier_tol", 1e-8)
            if verbose >= 2:
                ip_printer.print_problem_iter(state.niter,
                                              state.nfev,
                                              state.cg_niter,
                                              state.barrier_parameter,
                                              state.trust_radius,
                                              state.penalty,
                                              state.optimality,
                                              state.constr_violation)
            state.status = None
            if (callback is not None) and callback(state):
                state.status = 3
            elif state.optimality < gtol and state.constr_violation < gtol:
                state.status = 1
            elif state.trust_radius < xtol and state.barrier_parameter < barrier_tol:
                state.status = 2
            elif state.niter > max_iter:
                state.status = 0
            return state.status in (0, 1, 2, 3)

    if verbose >= 2:
        if method == 'equality_constrained_sqp':
            sqp_printer.print_header()
        if method == 'tr_interior_point':
            ip_printer.print_header()

    start_time = time.time()
    # Call inferior function to do the optimization
    if method == 'equality_constrained_sqp':
        result = equality_constrained_sqp(
            fun, grad, lagr_hess,
            constr.constr_eq,  constr.jac_eq,
            x0, fun0, grad0, constr_eq0, jac_eq0,
            stop_criteria, state, **options)
    elif method == 'tr_interior_point':
        result = tr_interior_point(
            fun, grad, lagr_hess,
            constr.n_ineq, constr.constr_ineq,
            constr.jac_ineq, constr.n_eq,
            constr.constr_eq, constr.jac_eq,
            x0, fun0, grad0, constr_ineq0, jac_ineq0,
            constr_eq0, jac_eq0, stop_criteria,
            constr.feasible_constr, sparse_jacobian,
            xtol, state, **options)

    result.execution_time = time.time() - start_time
    result.method = method
    result.message = TERMINATION_MESSAGES[result.status]

    if verbose >= 2:
        if method == 'equality_constrained_sqp':
            sqp_printer.print_footer()
        if method == 'tr_interior_point':
            ip_printer.print_footer()
    if verbose >= 1:
        print(result.message)
        print("Number of iteractions: {0}, function evaluations: {1}, "
              "CG iterations: {2}, optimality: {3:.2e}, "
              "constraint violation: {4:.2e}, execution time: {5:4.2} s."
              .format(result.niter, result.nfev, result.cg_niter,
                      result.optimality, result.constr_violation,
                      result.execution_time))
    return result
