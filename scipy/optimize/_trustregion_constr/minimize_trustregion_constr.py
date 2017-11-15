from __future__ import division, print_function, absolute_import
import numpy as np
from .._constraints import (NonlinearConstraint,
                            LinearConstraint,
                            BoxConstraint)
from .._canonical_constraint import (lagrangian_hessian,
                                    to_canonical,
                                    empty_canonical_constraint)
from .equality_constrained_sqp import equality_constrained_sqp
from .tr_interior_point import tr_interior_point
from warnings import warn
from copy import deepcopy
from scipy.sparse.linalg import LinearOperator
from .._quasi_newton_approx import QuasiNewtonApprox, BFGS, SR1
import scipy.sparse as sps
import time
from ..optimize import OptimizeResult
from .._numdiff import approx_derivative


TERMINATION_MESSAGES = {
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`xtol` termination condition is satisfied.",
    3: "`callback` function requested termination"
}


class ScalarFunction:
    """Define methods for evaluating a scalar function and its derivatives.

    This class define a scalar function F: R^n->R and methods for
    computing or approximating its first and second derivatives.
    """

    def __init__(self, fun, x0, args, grad='2-point', hess=BFGS(), finite_diff_options={}):
        self.x = np.atleast_1d(x0).astype(float)
        if grad in ('2-point', '3-point', 'cs'):
            finite_diff_options["method"] = grad
            self.x_diff = np.copy(self.x)
        if hess in ('2-point', '3-point', 'cs'):
            finite_diff_options["method"] = hess
            self.x_diff = np.copy(self.x)
        if grad in ('2-point', '3-point', 'cs') and \
           hess in ('2-point', '3-point', 'cs'):
            raise ValueError("Whenever the gradient is estimated via finite-differences, "
                             "we require the Hessian to be estimated using one of the "
                             "quasi-Newton strategies.")
        if isinstance(hess, QuasiNewtonApprox):
            self.x_prev = np.copy(x0)
            self.first_iteration = True

        # Define function
        self.f = fun(x0, *args)

        if grad in ('2-point', '3-point', 'cs'):
            def fun_wrapped(x):
                self.x_diff = x
                self.f = fun(x, *args)
                return self.f
        else:
            def fun_wrapped(x):
                return fun(x, *args)
        self.fun = fun_wrapped

        # Define gradient
        if callable(grad):
            self.g = np.atleast_1d(grad(x0, *args))
            if isinstance(hess, QuasiNewtonApprox):
                self.g_prev = np.copy(self.g)

            def grad_wrapped(x):
                return np.atleast_1d(grad(x, *args))

            if hess in ('2-point', '3-point', 'cs'):
                def grad_wrapped2(x):
                    self.x_diff = x
                    self.g = grad_wrapped(x)
                    return self.g

            elif isinstance(hess, QuasiNewtonApprox):
                def grad_wrapped2(x):
                    self.x_prev = self.x
                    self.g_prev = self.g
                    self.x = x
                    self.g = grad_wrapped(x)
                    return self.g
            else:
                grad_wrapped2 = grad_wrapped
        elif grad in ('2-point', '3-point', 'cs'):
            self.g = approx_derivative(fun, self.x, f0=self.f, **finite_diff_options)
            if isinstance(hess, QuasiNewtonApprox):
                self.g_prev = np.copy(self.g)

            def grad_wrapped(x):
                return approx_derivative(fun, x, f0=self.f, **finite_diff_options)

            if isinstance(hess, QuasiNewtonApprox):
                def grad_wrapped2(x):
                    if not np.array_equal(self.x_diff, x):
                        self.x_diff = x
                        self.f = fun(x, *args)
                    self.x_prev = self.x
                    self.g_prev = self.g
                    self.x = x
                    self.g = grad_wrapped(x)
                    return self.g

            else:
                def grad_wrapped2(x):
                    if not np.array_equal(self.x_diff, x):
                        self.x_diff = x
                        self.f = fun(x, *args)
                    return grad_wrapped(x)
        self.grad = grad_wrapped2

        # Define Hessian
        if callable(hess):
            self.H = hess(x0, *args)

            if sps.issparse(self.H):
                def hess_wrapped(x):
                    return sps.csr_matrix(hess(x, *args))
                self.H = sps.csr_matrix(self.H)

            elif isinstance(self.H, LinearOperator):
                def hess_wrapped(x):
                    return hess(x, *args)

            else:
                def hess_wrapped(x):
                    return np.atleast_2d(np.asarray(hess(x, *args)))
                self.H = np.atleast_2d(np.asarray(self.H))

        elif hess in ('2-point', '3-point', 'cs'):
            def hess_wrapped(x):
                if not np.array_equal(self.x_diff, x):
                    self.x_diff = x
                    self.g = grad_wrapped(x)
                return approx_derivative(grad_wrapped, x, f0=self.g, **finite_diff_options)
            self.H = approx_derivative(grad_wrapped, self.x, f0=self.g, **finite_diff_options)

        elif isinstance(hess, QuasiNewtonApprox):
            def hess_wrapped(x):
                if not np.array_equal(self.x, x):
                    self.x_prev = self.x
                    self.g_prev = self.g
                    self.x = x
                    self.g = grad_wrapped(x)
                delta_x = self.x - self.x_prev
                delta_grad = self.g - self.g_prev
                if self.first_iteration:
                    if np.linalg.norm(delta_x) != 0:
                        hess.instanciate_matrix(delta_x, delta_grad)
                        hess.scale_matrix(delta_x, delta_grad)
                        hess.update(delta_x, delta_grad)
                        self.first_iteration = False
                    else:
                        hess.instanciate_matrix(delta_x, delta_grad)
                else:
                    hess.update(delta_x, delta_grad)
                return hess
        else:
            hess_wrapped = None  
        self.hess = hess_wrapped


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
        print("|{0:>7}|{1:>7}|{2:>7}| {3:^1.2e} | {4:^1.2e} |"
              " {5:^1.2e} | {6:^1.2e} |"
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


def _minimize_trustregion_constr(fun, x0, args, grad,
                                 hess, constraints,
                                 xtol=1e-8, gtol=1e-8,
                                 barrier_tol=1e-8,
                                 sparse_jacobian=None,
                                 callback=None, maxiter=1000,
                                 verbose=0, finite_diff_options=None,
                                 initial_penalty=1.0, initial_trust_radius=1.0,
                                 initial_barrier_parameter=0.1,
                                 initial_tolerance=0.1,
                                 return_all=False, factorization_method=None,
                                 disp=False):
    """Minimize scalar function subject to constraints.

    Options
    -------
    xtol : float, optional
        Tolerance for termination by the change of the independent variable.
        The algorithm will terminate when ``delta < xtol``, where ``delta``
        is the algorithm trust-radius. Default is 1e-8.
    gtol : float, optional
        Tolerance for termination by the norm of the Lagrangian gradient.
        The algorithm will terminate when both the infinity norm (i.e. max
        abs value) of the Lagrangian gradient and the constraint violation
        are smaller than ``gtol``. Default is 1e-8.
    sparse_jacobian : {bool, None}
        The algorithm uses a sparse representation of the Jacobian if True
        and a dense representation if False. When sparse_jacobian is None
        the algorithm uses the more convenient option, using a sparse
        representation if at least one of the constraint Jacobians are sparse
        and a dense representation when they are all dense arrays.
    initial_trust_radius: float
        Initial trust-region radius. By defaut uses 1.0, as
        suggested in [1]_, p.19.
    initial_penalty : float
        Initial penalty for merit function. By defaut uses 1.0, as
        suggested in [1]_, p.19.
    initial_barrier_parameter: float
        Initial barrier parameter. Exclusive for 'tr_interior_point'
        method. By default uses 0.1, as suggested in [1]_, p. 19.
    initial_tolerance: float
        Initial subproblem tolerance. Exclusive for
        'tr_interior_point' method. By defaut uses 0.1,
        as suggested in [1]_, p. 19.
    return_all : bool, optional
        When True return the list of all vectors
        through the iterations.
    factorization_method : string, optional
        Method used for factorizing the jacobian matrix.
        Should be one of:

            - 'NormalEquation'.
            - 'AugmentedSystem'.
            - 'QRFactorization'.
            - 'SVDFactorization'.

        The factorization methods 'NormalEquation' and
        'AugmentedSystem' should be used only when
        ``sparse_jacobian=True``. The  projections
        required by the algorithm will be computed using,
        respectively, the the normal equation 
        and the augmented system approach explained in [1]_.
        'NormalEquation' computes the Cholesky
        factorization of ``(A A.T)`` and 'AugmentedSystem' 
        performes the LU factorization of an augmented system.
        They usually provide similar results.
        'NormalEquation' requires scikit-sparse
        installed. 'AugmentedSystem' is used by
        default for sparse matrices. 
        The methods 'QRFactorization'
        and 'SVDFactorization' should be used when
        ``sparse_jacobian=False``. They compute
        the required projections using, respectivelly,
        QR and SVD factorizations.
        The 'SVDFactorization' method can cope
        with Jacobian matrices with deficient row
        rank and will be used whenever other
        factorization methods fails (which may
        imply the conversion of sparse matrices
        to a dense format when required). By default uses
        'QRFactorization' for  dense matrices.
    finite_diff_options: dict, optional
        Dictionary with options to be passed on to `approx_derivative` method.
        These options will define the finite_difference approximation parameters.
    maxiter : int, optional
        Maximum number of algorithm iterations. By default ``maxiter=1000``
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

            * 0 (default) : work silently.
            * 1 : display a termination report.
            * 2 : display progress during iterations.
    disp : bool
        Set to True to force ``verbose`` to be greater or equal to 1,
        printing convergence mensages.

    Returns
    -------
    `OptimizeResult` with the following fields defined:
    x : ndarray, shape (n,)
        Solution found.
    s : ndarray, shape (n_ineq,)
        Slack variables at the solution. ``n_ineq`` is the total number
        of inequality constraints.
    v : ndarray, shape (n_ineq + n_eq,)
        Estimated Lagrange multipliers at the solution. ``n_ineq + n_eq``
        is the total number of equality and inequality constraints.
    niter : int
        Total number of iterations.
    nfev : int
        Total number of objective function evaluations.
    ngev : int
        Total number of objective function gradient evaluations.
    nhev : int
        Total number of Lagragian Hessian evaluations. Each time the
        Lagrangian Hessian is evaluated the objective function
        Hessian and the constraints Hessians are evaluated
        one time each.
    ncev : int
        Total number of constraint evaluations. The same couter
        is used for equality and inequality constraints, because
        they always are evaluated the same number of times.
    njev : int
        Total number of constraint Jacobian matrix evaluations.
        The same couter is used for equality and inequality
        constraint Jacobian matrices, because they always are
        evaluated the same number of times.
    cg_niter : int
        Total number of CG iterations.
    cg_info : Dict
        Dictionary containing information about the latest CG iteration:

            - 'niter' : Number of iterations.
            - 'stop_cond' : Reason for CG subproblem termination:

                1. Iteration limit was reached;
                2. Reached the trust-region boundary;
                3. Negative curvature detected;
                4. Tolerance was satisfied.

            - 'hits_boundary' : True if the proposed step is on the boundary
              of the trust region.

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
        Constraint violation at last iteration.
    optimality : float
        Norm of the Lagrangian gradient at last iteration.
    fun : float
        For the 'equality_constrained_sqp' method this is the objective
        function evaluated at the solution and for the 'tr_interior_point'
        method this is the barrier function evaluated at the solution.
    grad : ndarray, shape (n,)
        For the 'equality_constrained_sqp' method this is the gradient of the
        objective function evaluated at the solution and for the
        'tr_interior_point' method  this is the gradient of the barrier
        function evaluated at the solution.
    constr : ndarray, shape (n_ineq + n_eq,)
        For the 'equality_constrained_sqp' method this is the equality
        constraint evaluated at the solution and for the 'tr_interior_point'
        method this are the equality and inequality constraints evaluated at
        a given point (with the inequality constraints incremented by the
        value of the slack variables).
    jac : {ndarray, sparse matrix}, shape (n_ineq + n_eq, n)
        For the 'equality_constrained_sqp' method this is the Jacobian
        matrix of the equality constraint evaluated at the solution and
        for the tr_interior_point' method his is scaled augmented Jacobian
        matrix, defined as ``hat(A)`` in equation (19.36), reference [2]_,
        p. 581.
    """
    if disp and verbose == 0:
        verbose = 1
    if finite_diff_options is None:
        finite_diff_options = {}
    # Initial value
    if hess in ('2-point', '3-point', 'cs'):
        finite_diff_options["as_linear_operator"] = True
    objective = ScalarFunction(fun, x0, args, grad, hess, finite_diff_options)
    x0 = np.atleast_1d(x0).astype(float)
    n_vars = np.size(x0)

    # Put constraints in list format when needed
    if isinstance(constraints, (NonlinearConstraint,
                                LinearConstraint,
                                BoxConstraint)):
        constraints = [constraints]
    # Copy, evaluate and initialize constraints
    copied_constraints = [deepcopy(constr) for constr in constraints]
    for constr in copied_constraints:
        x0 = constr._evaluate_and_initialize(x0, sparse_jacobian)
    # Concatenate constraints
    if len(copied_constraints) == 0:
        constr = empty_canonical_constraint(x0, n_vars, sparse_jacobian)
    else:
        constr = to_canonical(copied_constraints)

    # Generate Lagrangian hess function
    lagr_hess = lagrangian_hessian(constr, objective.hess)

    # Construct OptimizeResult
    state = OptimizeResult(niter=0, nfev=1, ngev=1,
                           ncev=1, njev=1, nhev=0,
                           cg_niter=0, cg_info={})
    # Store values
    if return_all:
        state.allvecs = []
        state.allmult = []

    # Choose appropriate method
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
        def stop_criteria(state):
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

    if verbose >= 2:
        if method == 'equality_constrained_sqp':
            sqp_printer.print_header()
        if method == 'tr_interior_point':
            ip_printer.print_header()

    start_time = time.time()
    # Call inferior function to do the optimization
    if method == 'equality_constrained_sqp':
        if constr.n_ineq > 0:
            raise ValueError("'equality_constrained_sqp' does not "
                             "support inequality constraints.")

        def fun_and_constr(x):
            f = objective.fun(x)
            _, c_eq = constr.constr(x)
            return f, c_eq

        def grad_and_jac(x):
            g = objective.grad(x)
            _, J_eq = constr.jac(x)
            return g, J_eq

        result = equality_constrained_sqp(
            fun_and_constr, grad_and_jac, lagr_hess,
            x0, objective.f, objective.g,
            constr.c_eq0, constr.J_eq0,
            stop_criteria, state,
            initial_penalty, initial_trust_radius,
            return_all, factorization_method)

    elif method == 'tr_interior_point':
        if constr.n_ineq == 0:
            warn("The problem only has equality constraints. "
                 "The solver 'equality_constrained_sqp' is a "
                 "better choice for those situations.")
        result = tr_interior_point(
            objective.fun, objective.grad, lagr_hess,
            n_vars, constr.n_ineq, constr.n_eq,
            constr.constr, constr.jac,
            x0, objective.f, objective.g, constr.c_ineq0, constr.J_ineq0,
            constr.c_eq0, constr.J_eq0, stop_criteria,
            constr.enforce_feasibility,
            xtol, state, initial_barrier_parameter, initial_tolerance,
            initial_penalty, initial_trust_radius, return_all, factorization_method)
    else:
        raise ValueError("Unknown optimization ``method``.")

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
