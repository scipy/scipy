"""Trust-region interior point method.

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

from __future__ import division, print_function, absolute_import
import scipy.sparse as spc
import numpy as np
from .equality_constrained_sqp import equality_constrained_sqp
from scipy.sparse.linalg import LinearOperator

__all__ = ['tr_interior_point']


class BarrierSubproblem:
    """
    Barrier optimization problem:
        minimize fun(x) - barrier_parameter*sum(log(s))
        subject to: constr_eq(x)     = 0
                  constr_ineq(x) + s = 0
    """

    def __init__(self, x0, fun, grad, lagr_hess, n_vars, n_ineq, n_eq,
                 constr, jac, barrier_parameter, tolerance,
                 enforce_feasibility, global_stop_criteria,
                 xtol, fun0, grad0, constr_ineq0, jac_ineq0, constr_eq0,
                 jac_eq0):
        # Store parameters
        self.n_vars = n_vars
        self.x0 = x0
        self._fun = fun
        self._grad = grad
        self.lagr_hess = lagr_hess
        self._constr = constr
        self._jac = jac
        self.barrier_parameter = barrier_parameter
        self.tolerance = tolerance
        self.n_eq = n_eq
        self.n_ineq = n_ineq
        self.enforce_feasibility = enforce_feasibility
        self.global_stop_criteria = global_stop_criteria
        self.xtol = xtol
        # Auxiliar parameters
        self._x_f = x0
        self._f = fun0
        self._x_g = x0
        self._g = grad0
        self._x_c = x0
        self._c_ineq = constr_ineq0
        self._c_eq = constr_eq0
        self._x_j = x0
        self._J_ineq = jac_ineq0
        self._J_eq = jac_eq0

    def fun(self, x):
        if not np.array_equal(self._x_f, x):
            self._f = self._fun(x)
        return self._f

    def grad(self, x):
        if not np.array_equal(self._x_g, x):
            self._g = self._grad(x)
        return self._g

    def constr(self, x):
        if not np.array_equal(self._x_c, x):
            self._c_ineq, self._c_eq = self._constr(x)
        return self._c_ineq, self._c_eq

    def jac(self, x):
        if not np.array_equal(self._x_j, x):
            self._J_ineq, self._J_eq = self._jac(x)
        return self._J_ineq, self._J_eq

    def update(self, barrier_parameter, tolerance):
        self.barrier_parameter = barrier_parameter
        self.tolerance = tolerance

    def get_slack(self, z):
        return z[self.n_vars:self.n_vars+self.n_ineq]

    def get_variables(self, z):
        return z[:self.n_vars]

    def function(self, z):
        """Returns barrier function at given point.
        For z = [x, s], returns barrier function:
            function(z) = fun(x) - barrier_parameter*sum(log(s))
        """
        x = self.get_variables(z)
        s = self.get_slack(z)

        # Use technique from Nocedal and Wright book, ref [3]_, p.576,
        # to guarantee constraints from `enforce_feasibility`
        # stay feasible along iterations.
        c_ineq, _ = self.constr(x)
        s[self.enforce_feasibility] = -c_ineq[self.enforce_feasibility]
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
        c_ineq, c_eq = self.constr(x)
        return np.hstack((c_eq,
                          c_ineq + s))

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
        J_ineq, J_eq = self.jac(x)
        if self.n_ineq == 0:
            return J_eq
        else:
            if spc.issparse(J_eq) or spc.issparse(J_ineq):
                # It is expected that J_eq and J_ineq
                # are already `csr_matrix` because of
                # the way ``BoxConstraint``, ``NonlinearConstraint``
                # and ``LinearConstraint`` are defined.
                J_eq = spc.csr_matrix(J_eq)
                J_ineq = spc.csr_matrix(J_ineq)
                return self._assemble_sparse_jacobian(J_eq, J_ineq, s)
            else:
                S = np.diag(s)
                zeros = np.zeros((self.n_eq, self.n_ineq))
                # Convert to matrix
                if spc.issparse(J_ineq):
                    J_ineq = J_ineq.toarray()
                if spc.issparse(J_eq):
                    J_eq = J_eq.toarray()
                # Concatenate matrices
                return np.asarray(np.bmat([[J_eq, zeros],
                                           [J_ineq, S]]))

    def _assemble_sparse_jacobian(self, J_eq, J_ineq, s):
        """Assemble sparse jacobian given its components.

        Given ``J_eq``, ``J_ineq`` and ``s`` returns:
            jacobian = [ J_eq,     0     ]
                       [ J_ineq, diag(s) ]

        It is equivalent to:
            spc.bmat([[ J_eq,   None    ],
                      [ J_ineq, diag(s) ]], "csr")
        but significantly more efficient for this
        given structure.
        """
        n_vars, n_ineq, n_eq = self.n_vars, self.n_ineq, self.n_eq
        J_aux = spc.vstack([J_eq, J_ineq], "csr")
        indptr, indices, data = J_aux.indptr, J_aux.indices, J_aux.data
        new_indptr = indptr + np.hstack((np.zeros(n_eq, dtype=int),
                                         np.arange(n_ineq+1, dtype=int)))
        size = indices.size+n_ineq
        new_indices = np.empty(size)
        new_data = np.empty(size)
        mask = np.full(size, False, bool)
        mask[new_indptr[-n_ineq:]-1] = True
        new_indices[mask] = n_vars+np.arange(n_ineq)
        new_indices[~mask] = indices
        new_data[mask] = s
        new_data[~mask] = data
        J = spc.csr_matrix((new_data, new_indices, new_indptr),
                           (n_eq + n_ineq, n_vars + n_ineq))
        return J

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


def tr_interior_point(fun, grad, lagr_hess, n_vars, n_ineq, n_eq,
                      constr, jac, x0, fun0, grad0,
                      constr_ineq0, jac_ineq0, constr_eq0,
                      jac_eq0, stop_criteria,
                      enforce_feasibility, xtol, state,
                      initial_barrier_parameter=0.1,
                      initial_tolerance=0.1,
                      initial_penalty=1.0,
                      initial_trust_radius=1.0,
                      return_all=False,
                      factorization_method=None):
    """Trust-region interior points method.

    Solve problem:
        minimize fun(x)
        subject to: constr_ineq(x) <= 0
                    constr_eq(x) = 0
    using trust-region interior point method described in [1]_.
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

    # Default enforce_feasibility
    if enforce_feasibility is None:
        enforce_feasibility = np.zeros(n_ineq, bool)
    # Initial Values
    state.barrier_parameter = initial_barrier_parameter
    state.tolerance = initial_tolerance
    state.trust_radius = initial_trust_radius
    state.penalty = initial_penalty
    state.optimality = np.inf
    state.constr_violation = np.inf
    # Define barrier subproblem
    subprob = BarrierSubproblem(
        x0, fun, grad, lagr_hess, n_vars, n_ineq, n_eq, constr, jac,
        state.barrier_parameter, state.tolerance, enforce_feasibility,
        stop_criteria, xtol, fun0, grad0, constr_ineq0, jac_ineq0,
        constr_eq0, jac_eq0)

    # Define initial value for the slack variables
    s0 = np.maximum(-1.5*constr_ineq0, np.ones(n_ineq))
    # Define initial parameter for the first iteration.
    z = np.hstack((x0, s0))
    # Define trust region bounds
    trust_lb = np.hstack((np.full(subprob.n_vars, -np.inf),
                          np.full(subprob.n_ineq, -BOUNDARY_PARAMETER)))
    trust_ub = np.full(subprob.n_vars+subprob.n_ineq, np.inf)

    # If there are inequality constraints solve a
    # sequence of barrier problems
    first_barrier_prob = True
    while True:
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
        # Compute initial values
        fun0_subprob = subprob.function(z)
        grad0_subprob = subprob.gradient(z)
        constr0_subprob = subprob.constraints(z)
        jac0_subprob = subprob.jacobian(z)
        # Solve SQP subproblem
        state = equality_constrained_sqp(
            subprob.function, subprob.gradient,
            subprob.lagrangian_hessian, subprob.constraints,
            subprob.jacobian, z, fun0_subprob, grad0_subprob,
            constr0_subprob, jac0_subprob, subprob.stop_criteria,
            state, trust_lb, trust_ub, initial_penalty,
            state.trust_radius, subprob.scaling, return_all,
            factorization_method)
        z = state.x
        if stop_criteria(state):
            break

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
