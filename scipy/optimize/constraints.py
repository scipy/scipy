from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc


__all__ = ['CanonicalConstraint',
           'NonlinearConstraint',
           'LinearConstraint',
           'BoxConstraint',
           'parse_constraint',
           'concatenate_canonical_constraints']


class CanonicalConstraint:
    """Object representing canonical constraint

    Constraint of the form:
        constr_ineq(x) <= 0
        constr_eq(x) = 0

    Init Parameters
    ---------------
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
    hess : callable
        Hessian matrix:
            hess(x, v_eq, v_ineq) -> H

            - ``x``: array_like, shape (n,)
                Evaluation point.
            - ``v_eq``: array_like, shape (n_eq,), optional
                Lagrange multipliers for equality constraints.
            - ``v_ineq``: array_like, shape (n_ineq,), optional
                Lagrange multipliers for inequality constraints.
            - ``H``: LinearOperator (or sparse matrix or ndarray), shape (n, n)
                Lagrangian Hessian.
    """
    def __init__(self, n_ineq, constr_ineq, jac_ineq, n_eq, constr_eq, jac_eq, hess):
        self.n_ineq = n_ineq
        self.constr_ineq = constr_ineq
        self.jac_ineq = jac_ineq
        self.n_eq = n_eq
        self.constr_eq = constr_eq
        self.jac_eq = jac_eq
        self.hess = hess


class NonlinearConstraint:
    """Object representing nonlinear constraint

    Init Parameters
    ---------------
    fun : callable
        The function defining the constraint.
            fun(x) -> array_like, shape (m,)
    jac : callable
        Jacobian Matrix:
            jac(x) -> ndarray (or sparse matrix), shape (m, n)
        where x is a (n,) ndarray.
    hess : {callable, None}
        Hessian matrix of `dot(fun, v)`:
            hess(x, v) -> ndarray (or sparse matrix or LinearOperator), shape (n, n)
        where x is a (n,) ndarray and v is a (m,) ndarray.
    kind : tuple
        Specifies the type of contraint. Options for this
        parameters are:
            - ("interval", lb, ub): for a constraint of the type:
                lb <= fun(x) <= ub
            - ("greater", lb): for a constraint of the type:
                fun(x) >= lb
            - ("less", ub): for a constraint of the type:
                fun(x) <= ub
            - ("equals", c): for a constraint of the type:
                fun(x) == c
        where ``lb``,  ``ub`` and ``c`` are (m,) ndarrays and
        ``x`` is a (n,) ndarray.
    """
    def __init__(self, fun, jac, hess, kind):
        self.fun = fun
        self.jac = jac
        self.hess = hess
        self.kind = kind

    def to_canonical(self):
        # Parse constraints
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(self.kind)
        # Get dimensions
        n_eq = len(eq)
        n_ineq = len(ineq)
        # Set variables
        self.fun_x = None
        self.x = None
        self.jac_y = None
        self.y = None

        def constr_eq(x):
            if n_eq > 0:
                if not np.array_equal(x, self.x):
                    self.fun_x = self.fun(x)
                    self.x = x
                return self.fun_x[eq] - val_eq
            else:
                return np.empty((0,))

        def constr_ineq(x):
            if n_ineq > 0:
                if not np.array_equal(x, self.x):
                    self.fun_x = self.fun(x)
                    self.x = x
                return sign*(self.fun_x[ineq] - val_ineq)
            else:
                return np.empty((0,))

        def jac_eq(x):
            if n_eq > 0:
                if not np.array_equal(x, self.y):
                    self.jac_y = self.jac(x)
                    self.y = x
                return self.jac_y[eq, :]
            else:
                return np.empty((0, len(x)))

        def jac_ineq(x):
            if n_ineq > 0:
                if not np.array_equal(x, self.y):
                    self.jac_y = self.jac(x)
                    self.y = x
                if spc.issparse(self.jac_y):
                    D = spc.lil_matrix((n_ineq, n_ineq))
                    D.setdiag(sign)
                    return D*self.jac_y[ineq, :]
                else:
                    return np.multiply(self.jac_y[ineq, :],
                                       sign[:, np.newaxis])
            else:
                return np.empty((0, len(x)))

        if self.hess is None:
            hess = None
        else:
            def hess(x, v_eq=None, v_ineq=None):
                hess = self.hess
                v = np.zeros(fun_len)
                if v_eq is not None:
                    v[eq] += v_eq
                if v_ineq is not None:
                    v[ineq[sign == 1]] += v_ineq[sign == 1]
                    v[ineq[sign == -1]] -= v_ineq[sign == -1]
                return hess(x, v)

        return CanonicalConstraint(n_ineq, constr_ineq, jac_ineq,
                                   n_eq, constr_eq, jac_eq, hess)

class LinearConstraint:
    """Object representing linear constraint.

    Init Parameters
    ---------------
    A : ndarray (or sparse matrix), shape (m, n)
        Matrix for the linear constraint.
    kind : tuple
        Specifies the type of contraint. Options for this
        parameters are:
            - ("interval", lb, ub): for a constraint of the type:
                lb <= A x <= ub
            - ("greater", lb): for a constraint of the type:
                A x >= lb
            - ("less", ub): for a constraint of the type:
                A x <= ub
            - ("equals", c): for a constraint of the type:
                A x == c
        where ``lb``,  ``ub`` and ``c`` are (m,) ndarrays and
        ``x`` is a (n,) ndarray.
    """
    def __init__(self, A, kind):
        self.A = A
        self.kind = kind

    def to_nonlinear(self):
        def fun(x):
            return self.A.dot(x)

        def jac(x):
            return self.A

        return NonlinearConstraint(fun, jac, None, self.kind)

    def to_canonical(self):
        return self.to_nonlinear().to_canonical()

    


class BoxConstraint:
    """Object representing box constraint.

    Init Parameters
    ---------------
    kind : tuple
        Specifies the type of contraint. Options for this
        parameters are:
            - ("interval", lb, ub): for a constraint of the type:
                lb <= x <= ub
            - ("greater", lb): for a constraint of the type:
                x >= lb
            - ("less", ub): for a constraint of the type:
                x <= ub
        where ``lb``,  ``ub`` and ``c`` are (m,) ndarrays and
        ``x`` is a (n,) ndarray.
    """
    def __init__(self, kind):
        self.kind = kind

    def to_linear(self, sparse=True):
        _, _, _, _, _, fun_len = parse_constraint(self.kind)
        if sparse:
            A = spc.eye(fun_len).tocsc()
        else:
            A = np.eye(fun_len)
        return LinearConstraint(A, self.kind)

    def to_nonlinear(self, sparse=True):
        return self.to_linear(sparse).to_nonlinear()

    def to_canonical(self, sparse=True):
        return self.to_linear(sparse).to_canonical()


def parse_constraint(kind):
    """Read constraint type and return list of indices.

    Parameters
    ----------
    kind : tuple
        Specifies the type of contraint. Options for this
        parameters are:
            - ("interval", lb, ub): for a constraint of the type:
                lb[i] <= f[i] <= ub[i]
            - ("greater", lb): for a constraint of the type:
                f[i] >= lb[i]
            - ("less", ub): for a constraint of the type:
                f[i] <= ub[i]
            - ("equals", c): for a constraint of the type:
                f[i] == c[i] 
        where ``lb``,  ``ub`` and ``c`` are (m,) ndarrays.
    Returns
    -------
    eq : array_like
        A vector indicating equality constraints.
            len(eq) = number of equality constraints
    ineq : array_like
        A vector indicating inequality constraints.
            len(ineq) = number of inequality constraints
    val_eq : array_like
        Equality constraint right hand side:
            f[eq[i]] = val_eq[i]
    val_ineq : array_like
        Inequality constraint right hand side:
            sign[i]*(f[ineq[i]] - val_ineq[i]) <= 0
    sign : array_like
        Sign of inequality constraints.
    """

    if kind[0] == "equals":
        # Read values from input structure
        c = np.asarray(kind[1], dtype=float)
        # Set returns
        eq = np.arange(len(c), dtype=int)
        ineq = np.empty(0)
        val_eq = np.asarray(c)
        val_ineq = np.empty(0)
        sign = np.empty(0)
        fun_len = len(c)
    elif kind[0] in ("greater", "less", "interval"):
        # Constraint type
        if kind[0] == "greater":
            lb = np.asarray(kind[1], dtype=float)
            ub = np.full_like(lb, np.inf, dtype=float)
        elif kind[0] == "less":
            ub = np.asarray(kind[1], dtype=float)
            lb = np.full_like(ub, -np.inf, dtype=float)
        elif kind[0] == "interval":
            lb = np.asarray(kind[1], dtype=float)
            ub = np.asarray(kind[2], dtype=float)
            if len(lb) != len(ub):
                raise ValueError("Mismatching 'ub' and 'lb' dimensions.")
            if (lb > ub).any():
                raise ValueError("lb[i] > ub[i]")
        # Set auxiliar values
        arange = np.arange(len(lb), dtype=int)
        ones = np.ones(len(lb))
        lb_isinf = np.isinf(lb)
        ub_isinf = np.isinf(ub)
        eq_list = (lb == ub) & ~lb_isinf & ~ub_isinf
        # Set returns
        eq = arange[eq_list]
        val_eq = lb[eq_list]
        ineq = np.hstack((arange[~eq_list & ~lb_isinf],
                          arange[~eq_list & ~ub_isinf]))
        val_ineq = np.hstack((lb[~eq_list & ~lb_isinf],
                              ub[~eq_list & ~ub_isinf]))
        sign = np.hstack((-ones[~eq_list & ~lb_isinf],
                          ones[~eq_list & ~ub_isinf]))
        fun_len = len(lb)
    else:
        raise ValueError("kind '%s' not available." % kind[0])

    return eq, ineq, val_eq, val_ineq, sign, fun_len


def concatenate_canonical_constraints(constraints, sparse=True, hess=None):
    """Concatenate sequence of CanonicalConstraint's."""
    # Compute number of constraints
    n_eq = 0
    n_ineq = 0
    for constr in constraints:
        n_eq += constr.n_eq
        n_ineq += constr.n_ineq

    # Concatenate equality constraints
    def constr_eq(x):
        constr_eq_list = []
        for constr in constraints:
            constr_eq_list += [constr.constr_eq(x)]
        return np.hstack(constr_eq_list)

    # Concatenate inequality constraints
    def constr_ineq(x):
        constr_ineq_list = []
        for constr in constraints:
            constr_ineq_list += [constr.constr_ineq(x)]
        return np.hstack(constr_ineq_list)

    # Concatanate equality constraints Jacobian matrices
    def jac_eq(x):
        jac_eq_list = []
        for constr in constraints:
            J = constr.jac_eq(x)
            if not sparse and spc.issparse(J):
                jac_eq_list += [J.toarray()]
            else:
                jac_eq_list += [J]

        if sparse:
            return spc.vstack(jac_eq_list)
        else:
            return np.vstack(jac_eq_list)

    # Concatanate inequality constraints Jacobian matrices
    def jac_ineq(x):
        jac_ineq_list = []
        for constr in constraints:
            J = constr.jac_ineq(x)
            if not sparse and spc.issparse(J):
                jac_ineq_list += [J.toarray()]
            else:
                jac_ineq_list += [J]

        if sparse:
            return spc.vstack(jac_ineq_list)
        else:
            return np.vstack(jac_ineq_list)

    # Concatenate Hessians
    def lagr_hess(x, v_eq, v_ineq):
        n = len(x)
        if hess is not None:
            hess_list = [hess(x)]
        else:
            hess_list = []

        index_eq = 0
        index_ineq = 0
        for constr in constraints:
            if constr.hess is not None:
                hess_list += [constr.hess(x, v_eq[index_eq:index_eq+constr.n_eq],
                                          v_ineq[index_ineq:index_ineq+constr.n_ineq])]
            index_eq += constr.n_eq
            index_ineq += constr.n_ineq

        def matvec(p):
            result = np.zeros_like(p)
            for h in hess_list:
                result += h.dot(p)
            return result

        return spc.linalg.LinearOperator((n, n), matvec)


    return CanonicalConstraint(n_ineq, constr_ineq, jac_ineq,
                               n_eq, constr_eq, jac_eq, lagr_hess)
