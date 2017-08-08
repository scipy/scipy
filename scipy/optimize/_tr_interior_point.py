"""Trust-region interior points method"""

from __future__ import division, print_function, absolute_import
import scipy.sparse as spc
import numpy as np
from ._equality_constrained_sqp import equality_constrained_sqp
from .optimize import OptimizeResult
from .constraints import check_sparsity
from numpy.linalg import norm
from scipy.sparse.linalg import LinearOperator

__all__ = ['tr_interior_point']


class BarrierSubproblem:
    """
    Barrier optimization problem:
        minimize fun(x) - barrier_parameter*sum(log(s))
        subject to: constr_eq(x)     = 0
                  constr_ineq(x) + s = 0
    References
    ----------
    .. [1] Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal.
           "An interior point algorithm for large-scale nonlinear
           programming." SIAM Journal on Optimization 9.4 (1999): 877-900.
    .. [2] Byrd, Richard H., Guanghui Liu, and Jorge Nocedal.
           "On the local behavior of an interior point method for
           nonlinear programming." Numerical analysis 1997 (1997): 37-56.
    .. [3] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
           Second Edition (2006).
    """

    def __init__(self, x0, fun, grad, lagr_hess, n_ineq, constr_ineq,
                 jac_ineq, n_eq, constr_eq, jac_eq, barrier_parameter,
                 tolerance, feasible_constr_list, global_stop_criteria,
                 xtol, sparse_jacobian):
        # Compute number of variables
        self.n_vars, = np.shape(x0)
        # Store parameters
        self.x0 = x0
        self.fun = fun
        self.grad = grad
        self.lagr_hess = lagr_hess
        self._constr_ineq = constr_ineq
        self.jac_ineq = jac_ineq
        self._constr_eq = constr_eq
        self.jac_eq = jac_eq
        self.barrier_parameter = barrier_parameter
        self.tolerance = tolerance
        self.n_eq = n_eq
        self.n_ineq = n_ineq
        self.feasible_constr_list = feasible_constr_list
        self.global_stop_criteria = global_stop_criteria
        self.xtol = xtol
        self.sparse_jacobian = sparse_jacobian
        # Auxiliar parameters
        self._x_ineq = None
        self._x_eq = None
        self._c_ineq = None
        self._c_eq = None

    def constr_ineq(self, x):
        """Value of inequality constraint at current iteration.
        Avoid that multiple calls to constr_ineq cause
        the constraint to be evaluated multiple times."""
        if not np.array_equal(self._x_ineq, x):
            self._c_ineq = self._constr_ineq(x)
        return self._c_ineq

    def constr_eq(self, x):
        """Value of equality constraint at current iteration.
        Avoid that multiple calls to constr_ineq cause
        the constraint to be evaluated multiple times."""
        if not np.array_equal(self._x_eq, x):
            self._c_eq = self._constr_eq(x)
        return self._c_eq

    def update(self, barrier_parameter, tolerance):
        self.barrier_parameter = barrier_parameter
        self.tolerance = tolerance

    def get_slack(self, z):
        return z[self.n_vars:self.n_vars+self.n_ineq]

    def get_variables(self, z):
        return z[:self.n_vars]

    def s0(self):
        return np.ones(self.n_ineq)

    def z0(self):
        return np.hstack((self.x0, self.s0()))

    def function(self, z):
        """Returns barrier function at given point.
        For z = [x, s], returns barrier function:
            function(z) = fun(x) - barrier_parameter*sum(log(s))
        """
        x = self.get_variables(z)
        s = self.get_slack(z)

        # Use technique from Nocedal and Wright book, ref [3]_, p.576,
        # to guarantee constraints from `feasible_constr_list`
        # stay feasible along iterations.
        c_ineq = self.constr_ineq(x)
        s[self.feasible_constr_list] = -c_ineq[self.feasible_constr_list]
        log_s = [np.log(s_i) if s_i > 0 else -np.inf for s_i in s]

        return self.fun(x) - self.barrier_parameter*np.sum(log_s)

    def constraints(self, z):
        """Returns barrier problem constraints at given points.
        For z = [x, s], returns the constraints:
            constraints(z) = [   constr_eq(x)     ]
                             [ constr_ineq(x) + s ]
        """
        x = self.get_variables(z)
        s = self.get_slack(z)
        return np.hstack((self.constr_eq(x),
                          self.constr_ineq(x) + s))

    def scaling(self, z):
        """Returns scaling vector.
        Given by:
            scaling = [ones(n_vars), s]
        """
        s = self.get_slack(z)
        diag_elements = np.hstack((np.ones(self.n_vars), s))

        # Diagonal Matrix
        def matvec(vec):
            return diag_elements*vec
        return LinearOperator((self.n_vars+self.n_ineq,
                               self.n_vars+self.n_ineq),
                              matvec)

    def gradient(self, z):
        """Returns scaled gradient.
        Barrier scalled gradient  of the barrier problem
        by the previously defined scaling factor:
            gradient = [             grad(x)             ]
                       [ -barrier_parameter*ones(n_ineq) ]
        """
        x = self.get_variables(z)
        return np.hstack((self.grad(x),
                          -self.barrier_parameter*np.ones(self.n_ineq)))

    def jacobian(self, z):
        """Returns scaled Jacobian.
        Barrier scalled jacobian
        by the previously defined scaling factor:
            jacobian = [  jac_eq(x)  0  ]
                       [ jac_ineq(x) S  ]
        """
        x = self.get_variables(z)
        s = self.get_slack(z)
        jac_eq = self.jac_eq(x)
        jac_ineq = self.jac_ineq(x)
        if self.n_ineq == 0:
            return check_sparsity(jac_eq, self.sparse_jacobian)
        else:
            if self.sparse_jacobian is not None:
                sparse_jacobian = self.sparse_jacobian
            elif spc.issparse(jac_eq) or spc.issparse(jac_eq):
                sparse_jacobian = True
            else:
                sparse_jacobian = False

            if sparse_jacobian:
                S = spc.diags((s,), (0,))
                return spc.bmat([[jac_eq, None],
                                 [jac_ineq, S]], "csc")
            else:
                S = np.diag(s)
                zeros = np.zeros((self.n_eq, self.n_ineq))
                return np.asarray(np.bmat([[jac_eq, zeros],
                                           [jac_ineq, S]]))

    def lagrangian_hessian_x(self, z, v):
        """Returns Lagrangian Hessian (in relation to `x`) -> Hx"""
        x = self.get_variables(z)
        # Get lagrange multipliers relatated to nonlinear equality constraints
        v_eq = v[:self.n_eq]
        # Get lagrange multipliers relatated to nonlinear ineq. constraints
        v_ineq = v[self.n_eq:self.n_eq+self.n_ineq]
        lagr_hess = self.lagr_hess
        return lagr_hess(x, v_eq, v_ineq)

    def lagrangian_hessian_s(self, z, v):
        """Returns scaled Lagrangian Hessian (in relation to`s`) -> S Hs S"""
        s = self.get_slack(z)
        # Using the primal formulation:
        #     S Hs S = diag(s)*diag(barrier_parameter/s**2)*diag(s).
        # Reference [1]_ p. 882, formula (3.1)
        primal = self.barrier_parameter
        # Using the primal-dual formulation
        #     S Hs S = diag(s)*diag(v/s)*diag(s)
        # Reference [1]_ p. 883, formula (3.11)
        primal_dual = v[-self.n_ineq:]*s
        # Uses the primal-dual formulation for
        # positives values of v_ineq, and primal
        # formulation for the remaining ones.
        return np.where(v[-self.n_ineq:] > 0, primal_dual, primal)

    def lagrangian_hessian(self, z, v):
        """Returns scaled Lagrangian Hessian"""
        # Compute Hessian in relation to x and s
        Hx = self.lagrangian_hessian_x(z, v)
        if self.n_ineq > 0:
            S_Hs_S = self.lagrangian_hessian_s(z, v)

        # The scaled Lagragian Hessian is:
        #     [ Hx    0    ]
        #     [ 0   S Hs S ]
        def matvec(vec):
            vec_x = self.get_variables(vec)
            vec_s = self.get_slack(vec)
            if self.n_ineq > 0:
                return np.hstack((Hx.dot(vec_x), S_Hs_S*vec_s))
            else:
                return Hx.dot(vec_x)
        return LinearOperator((self.n_vars+self.n_ineq,
                               self.n_vars+self.n_ineq),
                              matvec)

    def stop_criteria(self, state):
        """Stop criteria to the barrier problem.
        The criteria here proposed is similar to formula (2.3)
        from [1]_, p.879.
        """
        return (state.optimality < self.tolerance
                and state.constr_violation < self.tolerance) \
            or self.global_stop_criteria(state) \
            or state.trust_radius < self.xtol


def tr_interior_point(fun, grad, lagr_hess, n_ineq, constr_ineq,
                      jac_ineq, n_eq, constr_eq, jac_eq, x0,
                      stop_criteria,
                      feasible_constr_list=None,
                      sparse_jacobian=None,
                      initial_barrier_parameter=0.1,
                      initial_tolerance=0.1,
                      initial_penalty=1.0,
                      initial_trust_radius=1.0,
                      return_all=False,
                      xtol_subproblem=1e-16):
    """Trust-region interior points method.

    Solve problem:
        minimize fun(x)
        subject to: constr_ineq(x) <= 0
                    constr_eq(x) = 0
    using trust-region interior point method described in [1]_.
    Parameters
    ----------
    fun : callable
        Objective function:
            fun(x) -> float
    grad : callable
        Gradient vector:
            grad(x) -> array_like, shape (n,)
    lagr_hess : callable
        Lagrangian hessian:
            hess(x, v_eq, v_ineq) -> H
            - ``x``: array_like, shape (n,)
                Evaluation point.
            - ``v_eq``: array_like, shape (n_eq,)
                Lagrange multipliers for equality constraints.
            - ``v_ineq``: array_like, shape (n_ineq,)
                Lagrange multipliers for inequality constraints.
            - ``H``: LinearOperator (or sparse matrix or ndarray), shape (n, n)
                Lagrangian Hessian.
    n_ineq : int
        Number of inequality constraints.
    constr_ineq : callable
        Inequality constraint:
            constr_ineq(x) -> array_like, shape (n_ineq,)
    jac_ineq : callable
        Inequality constraints Jacobian:
            jac_ineq(x) -> sparse matrix (or ndarray), shape (n_ineq, n)
    n_eq : int
        Number of equality constraints.
    constr_eq : callable
        Equality constraint:
            constr_eq(x) -> array_like, shape (n_eq,)
    jac_eq : callable
        Equality constraints Jacobian:
            jac_ineq(x) -> sparse matrix (or ndarray), shape (n_eq, n)
    x0 : array_like, shape (n,)
        Starting point.
    feasible_constr_list : array_like (boolean), shape (n_ineq,)
        List specifying inequality constraints. All the iterates generated
        by the optimization algorithm the algorithm will be feasible with
        respect to those constraints. It is important that the initial point
        ``x0`` respect the specified constraint, otherwise the algorithm will
        just fail.
    stop_criteria: callable
        Functions that returns True when stop criteria is fulfilled:
            stop_criteria(state)
    initial_tolerance: float
        Initial subproblem tolerance. By defaut uses 0.1.
    initial_barrier_parameter: float
        Initial barrier parameter. By defaut uses 0.1.
    initial_trust_radius: float
        Initial trust-region radius. By defaut uses 1.
    initial_penalty : float
        Initial penalty for merit function.
    return_all : bool, optional
        When ``true`` return the list of all vectors through the iterations.
    xtol : float, optional
        Tolerance for termination by the change of the independent variable
        for the barrier subproblem.
    Returns
    -------
    x : array_like, shape (n,)
        Solution to the optimization problem.
    state :
        Dictionary containing the following:
            - niter : Number of iterations.
            - trust_radius : Trust radius at last iteration.
            - v : Lagrange multipliers at the solution , shape (m,).
            - fun : Function evaluation at the solution.
            - grad : Gradient evaluation at the solution.
            - hess : Lagrangian Hessian at the solution.
            - constr : Constraints at the solution.
            - jac : Constraints jacobian at the solution.
            - opt : Optimality is the norm of gradient of the Lagrangian
              ``||grad L(x, v)||``, where ``grad L(x, v) = g(x) + A(x).T v``.
            - c_violation : Norm of the constraint violation ``||c(x)||``.
    References
    ----------
    .. [1] Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal.
           "An interior point algorithm for large-scale nonlinear
           programming." SIAM Journal on Optimization 9.4 (1999): 877-900.
    .. [2] Byrd, Richard H., Guanghui Liu, and Jorge Nocedal.
           "On the local behavior of an interior point method for
           nonlinear programming." Numerical analysis 1997 (1997): 37-56.
    .. [3] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
           Second Edition (2006).
    """
    # BOUNDARY_PARAMETER controls the decrease on the slack
    # variables. Represents ``tau`` from [1]_ p.885, formula (3.18).
    BOUNDARY_PARAMETER = 0.995
    # BARRIER_DECAY_RATIO controls the decay of the barrier parameter
    # and of the subproblem toloerance. Represents ``theta`` from [1]_ p.879.
    BARRIER_DECAY_RATIO = 0.2
    # TRUST_ENLARGEMENT controls the enlargement on trust radius
    # after each iteration
    TRUST_ENLARGEMENT = 5

    # Construct OptimizeResult
    state = OptimizeResult(niter=0, nfev=0, ngev=0,
                           ncev=0, njev=0, nhev=0,
                           cg_niter=0, cg_info={})
    # Default feasible_constr_list
    if feasible_constr_list is None:
        feasible_constr_list = np.zeros(n_ineq, bool)
    # Initial Values
    state.barrier_parameter = initial_barrier_parameter
    state.tolerance = initial_tolerance
    state.trust_radius = initial_trust_radius
    state.optimality = np.inf
    state.constr_violation = np.inf
    if return_all:
        state.allvecs = []
        state.allmult = []
    # Define barrier subproblem
    subprob = BarrierSubproblem(
        x0, fun, grad, lagr_hess, n_ineq, constr_ineq, jac_ineq,
        n_eq, constr_eq, jac_eq, state.barrier_parameter, state.tolerance,
        feasible_constr_list, stop_criteria, xtol_subproblem, sparse_jacobian)
    # Define initial parameter for the first iteration.
    z = subprob.z0()
    # Define trust region bounds
    trust_lb = np.hstack((np.full(subprob.n_vars, -np.inf),
                          np.full(subprob.n_ineq, -BOUNDARY_PARAMETER)))
    trust_ub = np.full(subprob.n_vars+subprob.n_ineq, np.inf)

    # If there are inequality constraints solve a
    # sequence of barrier problems
    first_barrier_prob = True
    while not stop_criteria(state):
        if not first_barrier_prob:
            # Update parameters
            state.trust_radius = max(initial_trust_radius,
                                     TRUST_ENLARGEMENT*state.trust_radius)
            # TODO: Use more advanced strategies from [2]_
            # to update this parameters.
            state.barrier_parameter *= BARRIER_DECAY_RATIO
            state.tolerance *= BARRIER_DECAY_RATIO
        first_barrier_prob = False
        # Update Barrier Problem
        subprob.update(state.barrier_parameter, state.tolerance)
        # Solve SQP subproblem
        state = equality_constrained_sqp(
            subprob.function,
            subprob.gradient,
            subprob.lagrangian_hessian,
            subprob.constraints,
            subprob.jacobian,
            z,
            subprob.stop_criteria,
            trust_lb,
            trust_ub,
            initial_penalty,
            state.trust_radius,
            subprob.scaling,
            state,
            return_all)
        z = state.x

    # Get x and s
    state.x = subprob.get_variables(z)
    state.s = subprob.get_slack(z)
    # Return all
    if return_all:
        allvecs = []
        allslack = []
        for z in state.allvecs:
            allvecs += [subprob.get_variables(z)]
            allslack += [subprob.get_slack(z)]
        state.allvecs = allvecs
        state.allslack = allslack

    return state
