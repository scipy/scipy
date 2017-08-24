from __future__ import division, print_function, absolute_import
import numpy as np
from ._constraints import (NonlinearConstraint,
                           LinearConstraint,
                           BoxConstraint)
from ._canonical_constraint import (lagrangian_hessian,
                                    to_canonical,
                                    empty_canonical_constraint)
from ._large_scale_constrained import (tr_interior_point,
                                       equality_constrained_sqp)
from warnings import warn
from copy import deepcopy
from scipy.sparse.linalg import LinearOperator
import scipy.sparse as sps
import time
from .optimize import OptimizeResult
from ._numdiff import approx_derivative


TERMINATION_MESSAGES = {
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`xtol` termination condition is satisfied.",
    3: "`callback` function requested termination"
}


class Memoize:
    "Memoize decorator, used to avoid repeated calls to the same function."
    def __init__(self, x0, f0):
        self._x = x0
        self._f = f0

    def __call__(self, function):
        def new_function(x):
            if not np.array_equal(x, self._x):
                self._x = x
                self._f = function(x)
            return self._f
        return new_function


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


def minimize_constrained(fun, x0, grad, hess='2-point', constraints=(),
                         method=None, xtol=1e-8, gtol=1e-8,
                         sparse_jacobian=None, options={},
                         callback=None, max_iter=1000,
                         verbose=0):
    """Minimize scalar function subject to constraints.

    Parameters
    ----------
    fun : callable
        The objective function to be minimized.

            ``fun(x) -> float``

        where x is an array with shape (n,).
    x0 : ndarray, shape (n,)
        Initial guess. Array of real elements of size (n,),
        where ``n`` is the number of independent variables.
    grad : callable
        Gradient of the objective function:

            ``grad(x) -> array_like, shape (n,)``

        where x is an array with shape (n,).
    hess : {callable, '2-point', '3-point', 'cs', None}, optional
        Method for computing the Hessian matrix. The keywords
        select a finite difference scheme for numerical
        estimation. The scheme '3-point' is more accurate, but requires
        twice as much operations compared to '2-point' (default). The
        scheme 'cs' uses complex steps, and while potentially the most
        accurate, it is applicable only when `fun` correctly handles
        complex inputs and can be analytically continued to the complex
        plane. If it is a callable, it should return the 
        Hessian matrix of `dot(fun, v)`:

            ``hess(x, v) -> {LinearOperator, sparse matrix, ndarray}, shape (n, n)``

        where x is a (n,) ndarray and v is a (m,) ndarray. When ``hess``
        is None it considers the hessian is a matrix filled with zeros.
    constraints : Constraint or List of Constraint's, optional
        A single object or a list of objects specifying
        constraints to the optimization problem.
        Available constraints are:

            - `BoxConstraint`
            - `LinearConstraint`
            - `NonlinearConstraint`

    method : {str, None}, optional
        Type of solver. Should be one of:

            - 'equality-constrained-sqp'
            - 'tr-interior-point'

        When None the more appropriate method is choosen.
        'equality-constrained-sqp' is chosen for problems that
        only have equality constraints and 'tr-interior-point'
        for general optimization problems.
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
    options : dict, optional
        A dictionary of solver options. Available options include:

            initial_trust_radius: float
                Initial trust-region radius. By defaut uses 1.0, as
                suggested in [1]_, p.19, immediatly after algorithm III.
            initial_penalty : float
                Initial penalty for merit function. By defaut uses 1.0, as
                suggested in [1]_, p.19, immediatly after algorithm III.
            initial_barrier_parameter: float
                Initial barrier parameter. Exclusive for 'tr_interior_point'
                method. By default uses 0.1, as suggested in [1]_ immediatly
                after algorithm III, p. 19.
            initial_tolerance: float
                Initial subproblem tolerance. Exclusive for
                'tr_interior_point' method. By defaut uses 0.1,
                as suggested in [1]_ immediatly after algorithm III, p. 19.
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
    callback : callable, optional
        Called after each iteration:

            ``callback(OptimizeResult state) -> bool``

        If callback returns True the algorithm execution is terminated.
        ``state`` is an `OptimizeResult` object, with the same fields
        as the ones from the return.
    max_iter : int, optional
        Maximum number of algorithm iterations. By default ``max_iter=1000``
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

            * 0 (default) : work silently.
            * 1 : display a termination report.
            * 2 : display progress during iterations.

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

    Notes
    -----
    Method 'equality_constrained_sqp' is an implementation of
    Byrd-Omojokun Trust-Region SQP method described [3]_ and
    in [2]_, p. 549. It solves equality constrained equality
    constrained optimization problems by solving, at each substep,
    a trust-region QP subproblem. The inexact solution of these
    QP problems using projected CG method makes this method
    appropriate for large-scale problems.

    Method 'tr_interior_point' is an implementation of the
    trust-region interior point method described in [1]_.
    It solves general nonlinear by introducing slack variables
    and solving a sequence of equality-constrained barrier problems
    for progressively smaller values of the barrier parameter.
    The previously described equality constrained SQP method is used
    to solve the subproblems with increasing levels of accuracy as
    the iterate gets closer to a solution. It is also an
    appropriate method for large-scale problems.

    References
    ----------
    .. [1] Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal.
           "An interior point algorithm for large-scale nonlinear
           programming." SIAM Journal on Optimization 9.4 (1999): 877-900.
    .. [2] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
           Second Edition (2006).
    .. [3] Lalee, Marucha, Jorge Nocedal, and Todd Plantega. "On the
           implementation of an algorithm for large-scale equality
           constrained optimization." SIAM Journal on
           Optimization 8.3 (1998): 682-706.
    """
    # Initial value
    x0 = np.atleast_1d(x0).astype(float)
    n_vars = np.size(x0)

    # Evaluate initial point
    f0 = fun(x0)
    g0 = np.atleast_1d(grad(x0))

    # Define Gradient
    if hess in ('2-point', '3-point', 'cs'):
        # Need to memoize gradient wrapper in order
        # to avoid repeated calls that occur
        # when using finite differences.
        @Memoize(x0, g0)
        def grad_wrapped(x):
            return np.atleast_1d(grad(x))

    else:
        def grad_wrapped(x):
            return np.atleast_1d(grad(x))

    # Check Hessian
    if callable(hess):
        H0 = hess(x0)

        if sps.issparse(H0):
            H0 = sps.csr_matrix(H0)

            def hess_wrapped(x):
                return sps.csr_matrix(hess(x))

        elif isinstance(H0, LinearOperator):
            def hess_wrapped(x):
                return hess(x)

        else:
            H0 = np.atleast_2d(np.asarray(H0))

            def hess_wrapped(x):
                return np.atleast_2d(np.asarray(hess(x)))

    elif hess in ('2-point', '3-point', 'cs'):
        approx_method = hess

        def hess_wrapped(x):
            return approx_derivative(grad_wrapped, x, approx_method,
                                     as_linear_operator=True)

    else:
        hess_wrapped = hess

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
    lagr_hess = lagrangian_hessian(constr, hess_wrapped)

    # Construct OptimizeResult
    state = OptimizeResult(niter=0, nfev=1, ngev=1,
                           ncev=1, njev=1, nhev=0,
                           cg_niter=0, cg_info={})
    # Store values
    return_all = options.get("return_all", False)
    if return_all:
        state.allvecs = []
        state.allmult = []

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
            elif (state.trust_radius < xtol
                  and state.barrier_parameter < barrier_tol):
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
        if constr.n_ineq > 0:
            raise ValueError("'equality_constrained_sqp' does not "
                             "support inequality constraints.")

        def fun_and_constr(x):
            f = fun(x)
            _, c_eq = constr.constr(x)
            return f, c_eq

        def grad_and_jac(x):
            g = grad_wrapped(x)
            _, J_eq = constr.jac(x)
            return g, J_eq

        result = equality_constrained_sqp(
            fun_and_constr, grad_and_jac, lagr_hess,
            x0, f0, g0, constr.c_eq0, constr.J_eq0,
            stop_criteria, state, **options)

    elif method == 'tr_interior_point':
        if constr.n_ineq == 0:
            warn("The problem only has equality constraints. "
                 "The solver 'equality_constrained_sqp' is a "
                 "better choice for those situations.")
        result = tr_interior_point(
            fun, grad_wrapped, lagr_hess,
            n_vars, constr.n_ineq, constr.n_eq,
            constr.constr, constr.jac,
            x0, f0, g0, constr.c_ineq0, constr.J_ineq0,
            constr.c_eq0, constr.J_eq0, stop_criteria,
            constr.enforce_feasibility,
            xtol, state, **options)
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
