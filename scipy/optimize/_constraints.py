from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import LinearOperator
from ._numdiff import approx_derivative
from warnings import warn
from ._quasi_newton_approx import QuasiNewtonApprox, BFGS, SR1


__all__ = ['NonlinearConstraint',
           'LinearConstraint',
           'BoxConstraint']


class VectorialFunction:
    """Define methods for evaluating a vectorial function and its derivatives.

    This class define a scalar function F: R^n->R^m and methods for
    computing or approximating its first and second derivatives.
    """

    def __init__(self, fun, x0, jac='2-point', hess=BFGS(),
                 finite_diff_options=None, sparse_jacobian=None):
        if finite_diff_options is None:
            finite_diff_options = {}
        self.x = np.atleast_1d(x0).astype(float)
        if jac in ('2-point', '3-point', 'cs'):
            finite_diff_options["method"] = jac
            self.x_diff = np.copy(self.x)
        if hess in ('2-point', '3-point', 'cs'):
            finite_diff_options["method"] = hess
            self.x_diff = np.copy(self.x)
        if jac in ('2-point', '3-point', 'cs') and \
           hess in ('2-point', '3-point', 'cs'):
            raise ValueError("Whenever the jacient is estimated via "
                             "finite-differences, we require the Hessian to "
                             "be estimated using one of the quasi-Newton "
                             "strategies.")
        if isinstance(hess, QuasiNewtonApprox):
            self.x_prev = np.copy(x0)
            self.first_iteration = True

        # Define function
        self.f = np.atleast_1d(fun(x0))

        def fun_wrapped(x):
                return np.atleast_1d(fun(x))

        if jac in ('2-point', '3-point', 'cs'):
            def fun_wrapped2(x):
                self.x_diff = x
                self.f = fun_wrapped(x)
                return self.f
        else:
            fun_wrapped2 = fun_wrapped
        self.fun = fun_wrapped2

        # Define jacobian
        if callable(jac):
            J0 = jac(x0)

            if sparse_jacobian or \
               (sparse_jacobian is None and sps.issparse(J0)):
                def jac_wrapped(x):
                    return sps.csr_matrix(jac(x))
                self.J = sps.csr_matrix(J0)
                if isinstance(hess, QuasiNewtonApprox):
                    self.J_prev = self.J.copy()
                self.sparse_jacobian = True

            elif sps.issparse(J0):
                def jac_wrapped(x):
                    return jac(x).toarray()
                self.J = J0.toarray()
                if isinstance(hess, QuasiNewtonApprox):
                    self.J_prev = np.copy(self.J)
                self.sparse_jacobian = False

            else:
                def jac_wrapped(x):
                    return np.atleast_2d(jac(x))
                self.J = np.atleast_2d(J0)
                if isinstance(hess, QuasiNewtonApprox):
                    self.J_prev = np.copy(self.J)
                self.sparse_jacobian = False

            if hess in ('2-point', '3-point', 'cs'):
                def jac_wrapped2(x):
                    self.x_diff = x
                    self.J = jac_wrapped(x)
                    return self.J

            elif isinstance(hess, QuasiNewtonApprox):
                def jac_wrapped2(x):
                    self.x_prev = self.x
                    self.J_prev = self.J
                    self.x = x
                    self.J = jac_wrapped(x)
                    return self.J

            else:
                jac_wrapped2 = jac_wrapped
        elif jac in ('2-point', '3-point', 'cs'):
            J0 = approx_derivative(fun, self.x, f0=self.f,
                                   **finite_diff_options)

            if sparse_jacobian or \
               (sparse_jacobian is None and sps.issparse(J0)):
                def jac_wrapped(x):
                    J = approx_derivative(fun, self.x, f0=self.f,
                                          **finite_diff_options)
                    return sps.csr_matrix(J)
                self.J = sps.csr_matrix(J0)
                if isinstance(hess, QuasiNewtonApprox):
                    self.J_prev = self.J.copy()
                self.sparse_jacobian = True

            elif sps.issparse(J0):
                def jac_wrapped(x):
                    J = approx_derivative(fun, self.x, f0=self.f,
                                          **finite_diff_options)
                    return J.toarray()
                self.J = J0.toarray()
                if isinstance(hess, QuasiNewtonApprox):
                    self.J_prev = np.copy(self.J)
                self.sparse_jacobian = False

            else:
                def jac_wrapped(x):
                    J = approx_derivative(fun, self.x, f0=self.f,
                                          **finite_diff_options)
                    return np.atleast_2d(J)
                self.J = np.atleast_2d(J0)
                if isinstance(hess, QuasiNewtonApprox):
                    self.J_prev = np.copy(self.J)
                self.sparse_jacobian = False

            if isinstance(hess, QuasiNewtonApprox):
                def jac_wrapped2(x):
                    if not np.array_equal(self.x_diff, x):
                        self.x_diff = x
                        self.f = fun_wrapped(x)
                    self.x_prev = self.x
                    self.J_prev = self.J
                    self.x = x
                    self.J = jac_wrapped(x)
                    return self.J

            else:
                def jac_wrapped2(x):
                    if not np.array_equal(self.x_diff, x):
                        self.x_diff = x
                        self.f = fun_wrapped(x)
                    return jac_wrapped(x)
        else:
            jac_wrapped = None
            jac_wrapped2 = None
        self.jac = jac_wrapped2

        # Define Hessian
        v0 = np.zeros_like(self.f)
        self.v_diff = v0
        if callable(hess):
            self.H = hess(x0, v0)

            if sps.issparse(self.H):
                def hess_wrapped(x, v):
                    return sps.csr_matrix(hess(x))
                self.H = sps.csr_matrix(self.H)

            elif isinstance(self.H, LinearOperator):
                def hess_wrapped(x, v):
                    return hess(x, v)

            else:
                def hess_wrapped(x, v):
                    return np.atleast_2d(np.asarray(hess(x, v)))
                self.H = np.atleast_2d(np.asarray(self.H))

        elif hess in ('2-point', '3-point', 'cs'):
            def jac_dot_v(x, v):
                J = jac_wrapped(x)
                return J.T.dot(v)
            self.H = approx_derivative(jac_dot_v, x0,
                                       f0=self.J.T.dot(v0), args=(v0,),
                                       **finite_diff_options)

            def hess_wrapped(x, v):
                if not np.array_equal(self.x_diff, x):
                    self.x_diff = x
                    self.J = jac_wrapped(x)
                return approx_derivative(jac_dot_v, x,
                                         f0=self.J.T.dot(v), args=(v,),
                                         **finite_diff_options)

        elif isinstance(hess, QuasiNewtonApprox):
            def hess_wrapped(x, v):
                if not np.array_equal(self.x, x):
                    self.x_prev = self.x
                    self.J_prev = self.J
                    self.x = x
                    self.J = jac_wrapped(x)
                delta_x = self.x - self.x_prev
                delta_grad = self.J.T.dot(v) - self.J_prev.T.dot(v)
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


class NonlinearConstraint:
    """Constraint imposing nonlinear relation on variables.

    Parameters
    ----------
    fun : callable
        The function defining the constraint.

            fun(x) -> array_like, shape (m,)

        where x is a (n,) ndarray and ``m``
        is the number of constraints.
    kind : {str, tuple}
        Specifies the type of contraint. To specify
        inequalities constraints use:

            - ``('greater',)`` for a constraint of the type:
                fun(x) >= 0
            - ``('less',)`` for a constraint of the type:
                fun(x) <= 0
            - ``('greater', lb)`` for a constraint of the type:
                fun(x) >= lb
            - ``('less', ub)`` for a constraint of the type:
                fun(x) <= ub

        To specify equalities constraints use:

            - ``('equals', c)`` for a constraint of the type:
                fun(x) == c
            - ``('equals',)`` for a constraint of the type:
                fun(x) == 0

        where ``lb``,  ``ub`` and ``c`` are (m,) ndarrays or
        scalar values. In the latter case, the same value
        will be repeated for all the constraints.
        Finally, bounds constraints can be specified using:

            - ``('interval', lb, ub)`` for a constraint of the type:
                lb <= fun(x) <= ub

        If ``lb[i] == ub[i]`` it define an equality constraint
        of the type ``fun(x)[i] == lb[i]``. Otherwise, the
        constraint requires the function to stay inside the
        specified boundaries. Use ``np.inf`` with an appropriate
        sign to disable inequalities constraints in one of the
        directions.
    jac : {callable,  '2-point', '3-point', 'cs'}, optional
        Method of computing the Jacobian matrix (an m-by-n matrix,
        where element (i, j) is the partial derivative of f[i] with
        respect to x[j]).  The keywords selects a finite difference
        scheme for numerical estimation of the jacobian matrix.
        Alternatively, if it is a callable, it should be a function that
        returns the matrix:

            ``jac(x) -> {ndarray, sparse matrix}, shape (m, n)``

        where x is a (n,) ndarray.
    hess : {callable, '2-point', '3-point', 'cs', None}
        Method for computing the Hessian matrix. The keywords
        select a finite difference scheme for numerical
        estimation.  Alternativelly,  a `QuasiNewtonApprox` object
        may be passed on, defining a quasi-Newton Hessian
        approximation method. If it is a callable, it should
        return the Hessian matrix of `dot(fun, v)`, that is:

            ``hess(x, v) -> {LinearOperator, sparse matrix, ndarray}, shape (n, n)``

        where x is a (n,) ndarray and v is a (m,) ndarray. During the solution
        of an optimization problem ``x`` and ``v`` will stand for,
        respectively, the vector of variables and the lagrange multipliers.
        When ``hess`` is None it considers the hessian is a matrix filled with zeros.
    enforce_feasibility : {list of bool, bool}, optional
        Specify wheather each of the constraint components needs to
        stay feasible. If ``True``  all the iterates generated by the
        optimization algorithm needs to be feasible in respect to the
        constraint, otherwise this will not be reinforced.
        A single value or a list are acceptable: A single value
        specifies the same behaviour for all constraints while
        a list would specify element-wise each constraints needs to
        stay feasible along the iterations and each does not.
        False by default.

    Notes
    -----
    Different methods are available for approximating the Hessian.
    The `QuasiNewtonApprox` object may define a quasi-Newton Hessian
    approximation. Available approximations are:

    - `BFGS`;
    - `SR1`;
    - `L-BFGS`.

    Finite difference schemes may be used for approximating
    either the gradient or the Hessian. We, however, do not allow
    its use for approximating both simultaneously. Hence
    whenever the gradient is estimated via finite-differences,
    we require the Hessian to be estimated using one of the
    quasi-Newton strategies.

    The scheme 'cs' is, potentially, the most accurate but
    it requires the function to correctly handles complex inputs
    and to be continuous in the complex plane. The scheme
    '3-point' is more accurate than '2-point' but requires twice
    as much operations.
    """
    def __init__(self, fun, kind, jac='2-point', hess=BFGS(),
                 enforce_feasibility=False,
                 finite_diff_options={}):
        self._fun = fun
        self.kind = kind
        self.finite_diff_options = finite_diff_options
        self._jac = jac
        self._hess = hess
        self.enforce_feasibility = enforce_feasibility
        self.isinitialized = False

    def _evaluate_and_initialize(self, x0, sparse_jacobian=None):
        if self._hess in ('2-point', '3-point', 'cs'):
            self.finite_diff_options["as_linear_operator"] = True
        constr = VectorialFunction(self._fun, x0, self._jac, self._hess,
                                   self.finite_diff_options, sparse_jacobian)

        self.fun = constr.fun
        self.jac = constr.jac
        self.hess = constr.hess
        self.sparse_jacobian = constr.sparse_jacobian
        self.x0 = constr.x
        self.f0 = constr.f
        self.J0 = constr.J
        self.n = constr.x.size
        self.m = constr.f.size
        self.kind = _check_kind(self.kind, self.m)
        self.enforce_feasibility = _check_enforce_feasibility(
            self.enforce_feasibility, self.m)
        if not _is_feasible(self.kind, self.enforce_feasibility, constr.f):
            raise ValueError("Unfeasible initial point. "
                             "Either set ``enforce_feasibility=False`` or "
                             "choose a new feasible initial point ``x0``.")

        self.isinitialized = True
        return constr.x


class LinearConstraint:
    """Linear constraint requiring ``A x`` to comply with a given relation.

    Parameters
    ----------
    A : {ndarray, sparse matrix}, shape (m, n)
        Matrix for the linear constraint.
    kind : {str, tuple}
        Specifies the type of contraint. To specify
        inequalities constraints use:

            - ``('greater',)`` for a constraint of the type:
                A x >= 0
            - ``('less',)`` for a constraint of the type:
                A x <= 0
            - ``('greater', lb)`` for a constraint of the type:
                A x >= lb
            - ``('less', ub)`` for a constraint of the type:
                A x <= ub

        To specify equalities constraints use:

            - ``('equals',)`` for a constraint of the type:
                A x == 0
            - ``('equals', c)`` for a constraint of the type:
                A x == c

        where ``lb``,  ``ub`` and ``c`` are (m,) ndarrays or
        scalar values. In the latter case, the same value
        will be repeated for all the constraints.
        Finally, bounds constraints can be specified using:

            - ``('interval', lb, ub)`` for a constraint of the type:
                lb <= A x <= ub

        If ``lb[i] == ub[i]`` it define an equality constraint
        of the type ``A x[i] == lb[i]``. Otherwise, the
        constraint requires the function to stay inside the
        specified boundaries. Use ``np.inf`` with an appropriate
        sign to disable inequalities constraints in one of the
        directions.
    enforce_feasibility : {list of bool, bool}, optional
        Specify wheather each of the constraint components needs to
        stay feasible. If ``True``  all the iterates generated by the
        optimization algorithm needs to be feasible in respect to the
        constraint, otherwise this will not be reinforced.
        A single value or a list are acceptable: A single value
        specifies the same behaviour for all constraints while
        a list would specify element-wise each constraints needs to 
        stay feasible along the iterations and each does not.
        False by default.
    """
    def __init__(self, A, kind, enforce_feasibility=False):
        self.A = A
        self.kind = kind
        self.enforce_feasibility = enforce_feasibility
        self.isinitialized = False

    def _evaluate_and_initialize(self, x0, sparse_jacobian=None):
        if sparse_jacobian or (sparse_jacobian is None
                               and sps.issparse(self.A)):
            self.A = sps.csr_matrix(self.A)
            self.sparse_jacobian = True
        elif sps.issparse(self.A):
            self.A = self.A.toarray()
            self.sparse_jacobian = False
        else:
            self.A = np.atleast_2d(self.A)
            self.sparse_jacobian = False

        x0 = np.atleast_1d(x0).astype(float)
        f0 = self.A.dot(x0)
        J0 = self.A

        self.x0 = x0
        self.f0 = f0
        self.J0 = J0
        self.n = x0.size
        self.m = f0.size
        self.kind = _check_kind(self.kind, self.m)
        self.enforce_feasibility = _check_enforce_feasibility(
            self.enforce_feasibility, self.m)
        if not _is_feasible(self.kind, self.enforce_feasibility, f0):
            raise ValueError("Unfeasible initial point. "
                             "Either set ``enforce_feasibility=False`` or "
                             "choose a new feasible initial point ``x0``.")

        self.isinitialized = True
        return x0

    def _to_nonlinear(self):
        if not self.isinitialized:
            raise RuntimeError("Trying to convert uninitialized constraint.")

        def fun(x):
            return self.A.dot(x)

        def jac(x):
            return self.A

        # Build Constraints
        nonlinear = NonlinearConstraint(fun, self.kind, jac, None,
                                        self.enforce_feasibility)
        nonlinear.isinitialized = True
        nonlinear.m = self.m
        nonlinear.n = self.n
        nonlinear.sparse_jacobian = self.sparse_jacobian
        nonlinear.fun = fun
        nonlinear.jac = jac
        nonlinear.hess = None
        nonlinear.x0 = self.x0
        nonlinear.f0 = self.f0
        nonlinear.J0 = self.J0
        return nonlinear


class BoxConstraint:
    """Box constraints confining varibles within a given interval.

    Parameters
    ----------
    kind : tuple
        Specifies the type of contraint. To specify
        inequalities constraints use:

            - ``('greater',)`` for a constraint of the type:
                x >= 0
            - ``('less',)`` for a constraint of the type:
                x <= 0
            - ``('greater', lb)`` for a constraint of the type:
                x >= lb
            - ``('less', ub)`` for a constraint of the type:
                x <= ub

        Bounds constraints can be specified using:

            - ``('interval', lb, ub)`` for a constraint of the type:
                lb <= x <= ub

        where ``lb`` and  ``ub`` are (m,) ndarrays or
        scalar values. In the latter case, the same value
        will be repeated for all the constraints. Use ``np.inf``
        with an appropriate sign to disable inequalities
        constraints in one of the  directions.
    enforce_feasibility : {list of bool, bool}, optional
        Specify wheather each of the constraint components needs to
        stay feasible. If ``True``  all the iterates generated by the
        optimization algorithm needs to be feasible in respect to the
        constraint, otherwise this will not be reinforced.
        A single value or a list are acceptable: A single value
        specifies the same behaviour for all constraints while
        a list would specify element-wise each constraints needs to 
        stay feasible along the iterations and each does not.
        False by default.
    """
    def __init__(self, kind, enforce_feasibility=False):
        self.kind = kind
        self.enforce_feasibility = enforce_feasibility
        self.isinitialized = False

    def _evaluate_and_initialize(self, x0, sparse_jacobian=None):
        x0 = np.atleast_1d(x0).astype(float)
        f0 = x0
        self.n = x0.size
        self.m = f0.size
        if sparse_jacobian or sparse_jacobian is None:
            J0 = sps.eye(self.n).tocsr()
            self.sparse_jacobian = True
        else:
            J0 = np.eye(self.n)
            self.sparse_jacobian = False

        self.J0 = J0
        self.kind = _check_kind(self.kind, self.m)
        self.enforce_feasibility = _check_enforce_feasibility(
            self.enforce_feasibility, self.m)
        self.isinitialized = True
        if not _is_feasible(self.kind, self.enforce_feasibility, f0):
            warn("The initial point was changed in order "
                 "to stay inside box constraints.")
            x0_new = _reinforce_box_constraint(self.kind,
                                               self.enforce_feasibility,
                                               x0)
            self.x0 = x0_new
            self.f0 = x0_new
            return x0_new
        else:
            self.x0 = x0
            self.f0 = f0
            return x0

    def _to_linear(self):
        if not self.isinitialized:
            raise RuntimeError("Trying to convert uninitialized constraint.")
        # Build Constraints
        linear = LinearConstraint(self.J0, self.kind,
                                  self.enforce_feasibility)
        linear.isinitialized = True
        linear.m = self.m
        linear.n = self.n
        linear.sparse_jacobian = self.sparse_jacobian
        linear.x0 = self.x0
        linear.f0 = self.f0
        linear.J0 = self.J0
        return linear

    def _to_nonlinear(self):
        if not self.isinitialized:
            raise RuntimeError("Trying to convert uninitialized constraint.")
        return self._to_linear()._to_nonlinear()


# ************************************************************ #
# **********           Auxiliary Functions           ********** #
# ************************************************************ #
def _check_kind(kind, m):
    if not isinstance(kind, (tuple, list, str)):
        raise ValueError("The parameter `kind` should be a tuple, "
                         " a list, or a string.")
    if isinstance(kind, str):
        kind = (kind,)
    if len(kind) == 0:
        raise ValueError("The parameter `kind` should not be empty.")

    n_args = len(kind)
    keyword = kind[0]
    if keyword not in ("greater", "less", "equals", "interval"):
        raise ValueError("Keyword `%s` not available." % keyword)
    if n_args in (1, 2) and keyword not in ("greater", "less", "equals") \
       or n_args == 3 and keyword not in ("interval"):
        raise ValueError("Invalid `kind` format.")
    if n_args == 1:
        kind = (keyword, 0)

    if keyword in ("greater", "less", "equals"):
        c = np.asarray(kind[1], dtype=float)
        if np.size(c) not in (1, m):
            if keyword == "greater":
                raise ValueError("`lb` has the wrong dimension.")
            if keyword == "less":
                raise ValueError("`ub` has the wrong dimension.")
            if keyword == "equals":
                raise ValueError("`c` has the wrong dimension.")
        c = np.resize(c, m)
        return (keyword, c)
    elif keyword == "interval":
        lb = np.asarray(kind[1], dtype=float)
        if np.size(lb) not in (1, m):
            raise ValueError("`lb` has the wrong dimension.")
        lb = np.resize(lb, m)
        ub = np.asarray(kind[2], dtype=float)
        if np.size(ub) not in (1, m):
            raise ValueError("`ub` has the wrong dimension.")
        ub = np.resize(ub, m)
        if (lb > ub).any():
            raise ValueError("lb[i] > ub[i].")
        return (keyword, lb, ub)


def _check_enforce_feasibility(enforce_feasibility, m):
    if isinstance(enforce_feasibility, bool):
        enforce_feasibility = np.full(m,
                                      enforce_feasibility,
                                      dtype=bool)
    else:
        enforce_feasibility = np.array(enforce_feasibility,
                                       dtype=bool)

        if enforce_feasibility.size != m:
            raise ValueError("The parameter 'enforce_feasibility' "
                             "has the wrong number of elements.")
    return enforce_feasibility


def _is_feasible(kind, enforce_feasibility, f0):
    keyword = kind[0]
    if keyword == "equals":
        lb = np.asarray(kind[1], dtype=float)
        ub = np.asarray(kind[1], dtype=float)
    elif keyword == "greater":
        lb = np.asarray(kind[1], dtype=float)
        ub = np.full_like(lb, np.inf, dtype=float)
    elif keyword == "less":
        ub = np.asarray(kind[1], dtype=float)
        lb = np.full_like(ub, -np.inf, dtype=float)
    elif keyword == "interval":
        lb = np.asarray(kind[1], dtype=float)
        ub = np.asarray(kind[2], dtype=float)
    else:
        raise RuntimeError("Never be here.")

    return ((lb[enforce_feasibility] <= f0[enforce_feasibility]).all()
            and (f0[enforce_feasibility] <= ub[enforce_feasibility]).all())


def _reinforce_box_constraint(kind, enforce_feasibility, x0,
                              relative_tolerance=0.01,
                              absolute_tolerance=0.01):
        """Move initial point ``x0`` to inside the box constraints."""
        x0 = np.copy(np.asarray(x0, dtype=float))
        keyword = kind[0]
        if keyword == "greater":
            lb = np.asarray(kind[1], dtype=float)
            ub = np.full_like(lb, np.inf, dtype=float)
        elif keyword == "less":
            ub = np.asarray(kind[1], dtype=float)
            lb = np.full_like(ub, -np.inf, dtype=float)
        elif keyword == "interval":
            lb = np.asarray(kind[1], dtype=float)
            ub = np.asarray(kind[2], dtype=float)

        x0_new = np.copy(x0)
        for i in range(np.size(x0)):
            if enforce_feasibility[i]:
                if not np.isinf(lb[i]):
                    lower_bound = min(lb[i]+absolute_tolerance,
                                      lb[i]+relative_tolerance*(ub[i]-lb[i]))
                    x0_new[i] = max(x0_new[i], lower_bound)
                if not np.isinf(ub[i]):
                    upper_bound = max(ub[i]-absolute_tolerance,
                                      ub[i]-relative_tolerance*(ub[i]-lb[i]))
                    x0_new[i] = min(x0_new[i], upper_bound)
        return x0_new
