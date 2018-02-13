from __future__ import division, print_function, absolute_import
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator
from .._differentiable_functions import VectorFunction
from .._constraints import (
    NonlinearConstraint, LinearConstraint, PreparedConstraint, strict_bounds)
from ..optimize import OptimizeResult
from .._differentiable_functions import ScalarFunction
from .equality_constrained_sqp import equality_constrained_sqp
from .canonical_constraint import (CanonicalConstraint,
                                   initial_constraints_as_canonical)
from .tr_interior_point import tr_interior_point
from .report import BasicReport, SQPReport, IPReport


TERMINATION_MESSAGES = {
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`xtol` termination condition is satisfied.",
    3: "`callback` function requested termination"
}


class HessianLinearOperator(object):
    """Build LinearOperator from hessp"""
    def __init__(self, hessp, n):
        self.hessp = hessp
        self.n = n

    def __call__(self, x, *args):
        def matvec(p):
            return self.hessp(x, p, *args)

        return LinearOperator((self.n, self.n), matvec=matvec)


class LagrangianHessian(object):
    """The Hessian of the Lagrangian as LinearOperator.

    The Lagrangian is computed as the objective function plus all the
    constraints multiplied with some numbers (Lagrange multipliers).
    """
    def __init__(self, n, objective_hess, constraints_hess):
        self.n = n
        self.objective_hess = objective_hess
        self.constraints_hess = constraints_hess

    def __call__(self, x, v_eq=np.empty(0), v_ineq=np.empty(0)):
        H_objective = self.objective_hess(x)
        H_constraints = self.constraints_hess(x, v_eq, v_ineq)

        def matvec(p):
            return H_objective.dot(p) + H_constraints.dot(p)

        return LinearOperator((self.n, self.n), matvec)


def update_state_sqp(state, objective, prepared_constraints, start_time,
                     trust_radius, penalty, cg_info):
    state.niter += 1
    state.x = objective.x
    state.fun = objective.f
    state.grad = objective.g
    state.nfev = objective.nfev
    state.njev = objective.ngev
    state.nhev = objective.nhev

    state.v = [c.fun.v for c in prepared_constraints]
    state.constr = [c.fun.f for c in prepared_constraints]
    state.jac = [c.fun.J for c in prepared_constraints]
    state.constr_nfev = [c.fun.nfev if isinstance(c.fun, VectorFunction) else 0
                         for c in prepared_constraints]
    state.constr_njev = [c.fun.njev if isinstance(c.fun, VectorFunction) else 0
                         for c in prepared_constraints]
    state.constr_nhev = [c.fun.nhev if isinstance(c.fun, VectorFunction) else 0
                         for c in prepared_constraints]

    state.execution_time = time.time() - start_time
    state.trust_radius = trust_radius
    state.penalty = penalty
    state.cg_niter += cg_info["niter"]
    state.cg_stop_cond = cg_info["stop_cond"]

    # Compute Lagrangian Gradient
    state.lagrangian_grad = np.copy(state.grad)
    for c in prepared_constraints:
        state.lagrangian_grad += c.fun.J.T.dot(c.fun.v)
    state.optimality = np.linalg.norm(state.lagrangian_grad, np.inf)

    # Compute maximum constraint violation
    state.constr_violation = 0
    for i in range(len(prepared_constraints)):
        lb = prepared_constraints[i].bounds[0]
        ub = prepared_constraints[i].bounds[1]
        c = state.constr[i]
        state.constr_violation = np.max([state.constr_violation,
                                         np.max(lb-c),
                                         np.max(c-ub)])

    return state


def update_state_ip(state, objective, prepared_constraints, start_time,
                    trust_radius, penalty, cg_info, barrier_parameter, tolerance):
    state = update_state_sqp(state, objective, prepared_constraints,
                             start_time, trust_radius, penalty, cg_info)
    state.barrier_parameter = barrier_parameter
    state.tolerance = tolerance
    return state


def _minimize_trustregion_constr(fun, x0, args, grad,
                                 hess, hessp, bounds, constraints,
                                 xtol=1e-8, gtol=1e-8,
                                 barrier_tol=1e-8,
                                 sparse_jacobian=None,
                                 callback=None, maxiter=1000,
                                 verbose=0, finite_diff_rel_step=None,
                                 initial_penalty=1.0, initial_trust_radius=1.0,
                                 initial_barrier_parameter=0.1,
                                 initial_tolerance=0.1,
                                 factorization_method=None,
                                 disp=False):
    """Minimize a scalar function subject to constraints.

    Options
    -------
    gtol : float, optional
        Tolerance for termination by the norm of the Lagrangian gradient.
        The algorithm will terminate when both the infinity norm (i.e. max
        abs value) of the Lagrangian gradient and the constraint violation
        are smaller than ``gtol``. Default is 1e-8.
    xtol : float, optional
        Tolerance for termination by the change of the independent variable.
        The algorithm will terminate when ``tr_radius < xtol``, where
        ``tr_radius`` is the radius of the trust region used in the algorithm.
        Default is 1e-8.
    barrier_tol : float, optional
        Threshold on the barrier parameter for the algorithm termination.
        When inequality constraints are present the algorithm will terminated
        only when the barrier parameter is less than `barrier_tol`.
        Default is 1e-8.
    sparse_jacobian : {bool, None}, optional
        Determines how to represent the Jacobian of the constraints. If bool,
        then Jacobians of all constraints will be converted to the
        corresponding format. If None (default), then Jacobians won't be
        converted, but the algorithm can proceed only if they all have the
        same format.
    initial_trust_radius: float, optional
        Initial radius of the trust region.
        Default is 1 (recommended in [1]_, p. 19).
    initial_penalty : float, optional
        Initial penalty for the merit function.
        Default is 1 (recommended in [1]_, p 19).
    initial_barrier_parameter: float, optional
        Initial barrier parameter. Used only when inequality constraints
        are present. Default is 0.1 (recommended in [1]_ p. 19).
    initial_tolerance: float, optional
        Initial subproblem tolerance. Used only when inequality constraints
        are present. Default is 0.1 (recommended in [1]_ p. 19).
    factorization_method : string or None, optional
        Method to factorize the Jacobian of the constraints.
        Use None (default) for the auto selection or one of:

            - 'NormalEquation'
            - 'AugmentedSystem'
            - 'QRFactorization'
            - 'SVDFactorization'

        The factorization methods 'NormalEquation' and 'AugmentedSystem'
        should be used only when ``sparse_jacobian=True``. The projections
        required by the algorithm will be computed using, respectively,
        the the normal equation and the augmented system approach explained
        in [1]_. 'NormalEquation' computes the Cholesky factorization of
        ``(A A.T)`` and 'AugmentedSystem' performs the LU factorization
        of an augmented system. They usually provide similar results.
        'NormalEquation' requires scikit-sparse installed. 'AugmentedSystem'
        is used by default for sparse matrices. The methods 'QRFactorization'
        and 'SVDFactorization' should be used when ``sparse_jacobian=False``.
        They compute the required projections using, respectively,
        QR and SVD factorizations. The 'SVDFactorization' method can cope
        with Jacobian matrices with deficient row rank and will be used
        whenever other factorization methods fail (which may imply the
        conversion of sparse matrices to a dense format when required).
        By default uses 'QRFactorization' for  dense matrices.
    finite_diff_rel_step : None or array_like, optional
        Relative step size for the finite difference approximation.
    maxiter : int, optional
        Maximum number of algorithm iterations. Default is 1000.
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

            * 0 (default) : work silently.
            * 1 : display a termination report.
            * 2 : display progress during iterations.
            * 3 : display progress during iterations (more complete report).

    disp : bool, optional
        If True (default) then `verbose` will be set to 1 if it was 0.

    Returns
    -------
    `OptimizeResult` with the following fields defined:
    x : ndarray, shape (n,)
        Solution found.
    fun : float
        Objective function at the solution.
    grad : ndarray, shape (n,)
        Gradient of the objective function at the solution.
    lagrangian_grad : ndarray, shape (n,)
        Gradient of the Lagrangian function at the solution.
    nfev : integer
        Number of objective function evaluations.
    ngev : integer
        Number of gradient evaluations. Will be zero when
        using finite-differences approximations.
    nhev : integer
        Number of Hessian evaluations. Will be zero when
        using finite-differences or quasi-Newton approximations.
    niter : int
        Total number of iterations.
    cg_niter : int
        Total number of CG iterations.
    cg_stop_cond : int
        Reason for CG subproblem termination at last iteration:

            * 0 : CG subproblem not evaluated;
            * 1 : Iteration limit was reached;
            * 2 : Reached the trust-region boundary;
            * 3 : Negative curvature detected;
            * 4 : Tolerance was satisfied.

    execution_time : float
        Total execution time.
    trust_radius : float
        Trust radius at the last iteration.
    penalty : float
        Penalty function at last iteration.
    tolerance : float
        Tolerance for barrier subproblem at the last iteration.
        Exclusive for 'tr_interior_point'.
    barrier_parameter : float
        Barrier parameter at the last iteration. Exclusive for
        'tr_interior_point'.
    status : {0, 1, 2, 3}
        Termination status:

            * 0 : The maximum number of function evaluations is exceeded.
            * 1 : `gtol` termination condition is satisfied.
            * 2 : `xtol` termination condition is satisfied.
            * 3 : `callback` function requested termination.

    message : str
        Termination message.
    method : {'equality_constrained_sqp', 'tr_interior_point'}
        Optimization method used.
    constr_violation : float
        Maximum constraint violation at the solution.
    optimality : float
        Infinity norm of the Lagrangian gradient at the solution.
    constr : list of ndarray
        List of constraint values at the solution. This list will give
        the values in the same order used  in `constraints`. When bound
        constraints are present one extra element will be include(in
        the end of the list) to account to the values corresponding
        to these constraints
    jac : list of {ndarray, sparse matrix}
        List of constraints Jacobian matrices evaluated at the solution.
        This list will give the values in the same order used in
        `constraints`. When bound constraints are present one extra element
        will be include (in the end of the list) to account to the values
        of these constraints Jacobian matrices.
    v : list of ndarray
        List of estimated Lagrange multipliers at the solution.
        This list will give the lagrange multiplers in the same
        order  used in `constraints`. When bound constraints are
        present one extra element will be included (in the end of
        the list) to account to the Lagrange multipliers of these
        constraints.
    constr_nfev : list of int
        Number of constraint evaluations for each of the constraints.
        This list will give the values in the same order used in
        `constraints`. When bound constraints are present one extra element
        will be include (in the end of the list) to account to the values
        of these constraints Jacobian matrices.
    constr_njev : list of int
        Number of Jacobian matrix evaluations for each of the constraints.
        This list will give the values in the same order used in
        `constraints`. When bound constraints are present one extra element
        will be include (in the end of the list) to account to the values
        of these constraints Jacobian matrices. Counters associated
        with linear and bound constraints or with constraints using
        finite-difference aproximations will be zero.
    constr_nhev : list of int
        Number of Hessian evaluations for each of the constraints.
        This list will give the values in the same order used in
        `constraints`. When bound constraints are present one extra element
        will be include (in the end of the list) to account to the values
        of these constraints Jacobian matrices. Counters associated
        with linear and bound constraints or with constraints using
        finite-difference or quasi-Newton aproximations will be zero.
    """
    x0 = np.atleast_1d(x0).astype(float)
    n_vars = np.size(x0)
    if callable(hessp) and hess is None:
        hess = HessianLinearOperator(hessp, n_vars)
    if disp and verbose == 0:
        verbose = 1

    if bounds is not None:
        finite_diff_bounds = strict_bounds(bounds.lb, bounds.ub,
                                           bounds.keep_feasible, n_vars)
    else:
        finite_diff_bounds = (-np.inf, np.inf)

    # Define Objective Funciton
    objective = ScalarFunction(fun, x0, args, grad, hess,
                               finite_diff_rel_step, finite_diff_bounds)

    # Put constraints in list format when needed
    if isinstance(constraints, (NonlinearConstraint, LinearConstraint)):
        constraints = [constraints]

    # Prepare constraints.
    prepared_constraints = [
        PreparedConstraint(c, x0, sparse_jacobian, finite_diff_bounds)
        for c in constraints]

    # Check that all constraints are either sparse or dense.
    n_sparse = sum(c.fun.sparse_jacobian for c in prepared_constraints)
    if 0 < n_sparse < len(prepared_constraints):
        raise ValueError("All constraints must have the same kind of the "
                         "Jacobian --- either all sparse or all dense. "
                         "You can set the sparsity globally by setting "
                         "`sparse_jacobian` to either True of False.")
    if prepared_constraints:
        sparse_jacobian = n_sparse > 0

    if bounds is not None:
        prepared_constraints.append(PreparedConstraint(bounds, x0,
                                                       sparse_jacobian))

    # Concatenate initial constraints to the canonical form.
    c_eq0, c_ineq0, J_eq0, J_ineq0 = initial_constraints_as_canonical(
        n_vars, prepared_constraints, sparse_jacobian)

    # Prepare all canonical constraints and concatenate it into one.
    canonical_all = [CanonicalConstraint.from_PreparedConstraint(c)
                     for c in prepared_constraints]

    if len(canonical_all) == 0:
        canonical = CanonicalConstraint.empty(n_vars)
    elif len(canonical_all) == 1:
        canonical = canonical_all[0]
    else:
        canonical = CanonicalConstraint.concatenate(canonical_all,
                                                    sparse_jacobian)

    # Generate the Hessian of the Lagrangian.
    lagrangian_hess = LagrangianHessian(n_vars, objective.hess, canonical.hess)

    # Choose appropriate method
    if canonical.n_ineq == 0:
        method = 'equality_constrained_sqp'
    else:
        method = 'tr_interior_point'

    # Construct OptimizeResult
    state = OptimizeResult(
        niter=0, nfev=0, njev=0, nhev=0,
        cg_niter=0, cg_stop_cond=0,
        fun=objective.f, grad=objective.g,
        lagrangian_grad=np.copy(objective.g),
        constr=[c.fun.f for c in prepared_constraints],
        jac=[c.fun.J for c in prepared_constraints],
        constr_nfev=[0 for c in prepared_constraints],
        constr_njev=[0 for c in prepared_constraints],
        constr_nhev=[0 for c in prepared_constraints],
        v=[c.fun.v for c in prepared_constraints],
        method=method)

    # Start counting
    start_time = time.time()

    # Define stop criteria
    if method == 'equality_constrained_sqp':
        def stop_criteria(state, optimality, constr_violation, trust_radius,
                          penalty, cg_info):
            state = update_state_sqp(state, objective, prepared_constraints,
                                     start_time, trust_radius, penalty, cg_info)
            if verbose == 2:
                BasicReport.print_iteration(state.niter,
                                            state.nfev,
                                            state.cg_niter,
                                            state.fun,
                                            state.trust_radius,
                                            state.optimality,
                                            state.constr_violation)
            elif verbose > 2:
                SQPReport.print_iteration(state.niter,
                                          state.nfev,
                                          state.cg_niter,
                                          state.fun,
                                          state.trust_radius,
                                          state.optimality,
                                          state.constr_violation,
                                          state.penalty,
                                          state.cg_stop_cond)
            state.status = None
            if (callback is not None) and callback(state.x, state):
                state.status = 3
            elif state.optimality < gtol and state.constr_violation < gtol:
                state.status = 1
            elif state.trust_radius < xtol:
                state.status = 2
            elif state.niter > maxiter:
                state.status = 0
            return state.status in (0, 1, 2, 3)
    elif method == 'tr_interior_point':
        def stop_criteria(state, trust_radius, penalty, cg_info, barrier_parameter,
                          tolerance):
            state = update_state_ip(state, objective, prepared_constraints,
                                    start_time, trust_radius, penalty, cg_info,
                                    barrier_parameter, tolerance)
            if verbose == 2:
                BasicReport.print_iteration(state.niter,
                                            state.nfev,
                                            state.cg_niter,
                                            state.fun,
                                            state.trust_radius,
                                            state.optimality,
                                            state.constr_violation)
            elif verbose > 2:
                IPReport.print_iteration(state.niter,
                                         state.nfev,
                                         state.cg_niter,
                                         state.fun,
                                         state.trust_radius,
                                         state.optimality,
                                         state.constr_violation,
                                         state.penalty,
                                         state.barrier_parameter,
                                         state.cg_stop_cond)
            state.status = None
            if (callback is not None) and callback(state.x, state):
                state.status = 3
            elif state.optimality < gtol and state.constr_violation < gtol:
                state.status = 1
            elif (state.trust_radius < xtol
                  and state.barrier_parameter < barrier_tol):
                state.status = 2
            elif state.niter > maxiter:
                state.status = 0
            return state.status in (0, 1, 2, 3)

    if verbose == 2:
        BasicReport.print_header()
    elif verbose > 2:
        if method == 'equality_constrained_sqp':
            SQPReport.print_header()
        elif method == 'tr_interior_point':
            IPReport.print_header()

    # Call inferior function to do the optimization
    if method == 'equality_constrained_sqp':
        if canonical.n_ineq > 0:
            raise ValueError("'equality_constrained_sqp' does not "
                             "support inequality constraints.")

        def fun_and_constr(x):
            f = objective.fun(x)
            c_eq, _ = canonical.fun(x)
            return f, c_eq

        def grad_and_jac(x):
            g = objective.grad(x)
            J_eq, _ = canonical.jac(x)
            return g, J_eq

        _, result = equality_constrained_sqp(
            fun_and_constr, grad_and_jac, lagrangian_hess,
            x0, objective.f, objective.g,
            c_eq0, J_eq0,
            stop_criteria, state,
            initial_penalty, initial_trust_radius,
            factorization_method)

    elif method == 'tr_interior_point':
        _, result = tr_interior_point(
            objective.fun, objective.grad, lagrangian_hess,
            n_vars, canonical.n_ineq, canonical.n_eq,
            canonical.fun, canonical.jac,
            x0, objective.f, objective.g,
            c_ineq0, J_ineq0, c_eq0, J_eq0,
            stop_criteria,
            canonical.keep_feasible,
            xtol, state, initial_barrier_parameter, initial_tolerance,
            initial_penalty, initial_trust_radius,
            factorization_method)
    else:
        raise ValueError("Unknown optimization ``method``.")

    result.message = TERMINATION_MESSAGES[result.status]

    if verbose == 2:
        BasicReport.print_footer()
    elif verbose > 2:
        if method == 'equality_constrained_sqp':
            SQPReport.print_footer()
        if method == 'tr_interior_point':
            IPReport.print_footer()
    if verbose >= 1:
        print(result.message)
        print("Number of iterations: {}, function evaluations: {}, "
              "CG iterations: {}, optimality: {:.2e}, "
              "constraint violation: {:.2e}, execution time: {:4.2} s."
              .format(result.niter, result.nfev, result.cg_niter,
                      result.optimality, result.constr_violation,
                      result.execution_time))
    return result
