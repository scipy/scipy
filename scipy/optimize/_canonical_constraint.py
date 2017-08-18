from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc
from ._constraints import (check_sparsity)

__all__ = ['CanonicalConstraint',
           'nonlinear_to_canonical',
           'linear_to_canonical',
           'box_to_canonical',
           'parse_constraint',
           'empty_canonical_constraint',
           'generate_lagrangian_hessian',
           'concatenate_canonical_constraints']


class CanonicalConstraint:
    """Canonical constraint

    Constraint of the form:
        constr_ineq(x) <= 0
        constr_eq(x) = 0

    Parameters
    ----------
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
            - ``H``: {LinearOperator, sparse matrix, ndarray}, shape (n, n)
                Hessian matrix.

    enforce_feasibility : list of boolean, shape (n_ineq,)
        List of booleans containing ``True`` if the correspondent inequality
        constraint must be feasible along the iteractions and ``False``
        otherwise.
    """
    def __init__(self, n_ineq, constr_ineq, jac_ineq,
                 n_eq, constr_eq, jac_eq, hess, enforce_feasibility):
        self.n_ineq = n_ineq
        self.constr_ineq = constr_ineq
        self.jac_ineq = jac_ineq
        self.n_eq = n_eq
        self.constr_eq = constr_eq
        self.jac_eq = jac_eq
        self.hess = hess
        self.enforce_feasibility = enforce_feasibility


def nonlinear_to_canonical(nonlinear, sparse_jacobian=None):
    # Parse constraints
    eq, ineq, val_eq, val_ineq, sign, fun_len \
        = parse_constraint(nonlinear.kind)
    # Get dimensions
    n_eq = len(eq)
    n_ineq = len(ineq)
    # Set variables
    nonlinear.fun_x = None
    nonlinear.x = None
    nonlinear.jac_y = None
    nonlinear.y = None

    def constr_eq(x):
        if n_eq > 0:
            if not np.array_equal(x, nonlinear.x):
                nonlinear.fun_x = nonlinear.fun(x)
                nonlinear.x = x
            return nonlinear.fun_x[eq] - val_eq
        else:
            return np.empty((0,))

    def constr_ineq(x):
        if n_ineq > 0:
            if not np.array_equal(x, nonlinear.x):
                nonlinear.fun_x = nonlinear.fun(x)
                nonlinear.x = x
            return sign*(nonlinear.fun_x[ineq] - val_ineq)
        else:
            return np.empty((0,))

    def jac_eq(x):
        if n_eq > 0:
            if not np.array_equal(x, nonlinear.y):
                nonlinear.jac_y = nonlinear.jac(x)
                nonlinear.y = x
            J = nonlinear.jac_y[eq, :]
            return check_sparsity(J, sparse_jacobian)
        else:
            return np.empty((0, len(x)))

    def jac_ineq(x):
        if n_ineq > 0:
            if not np.array_equal(x, nonlinear.y):
                nonlinear.jac_y = nonlinear.jac(x)
                nonlinear.y = x
            if spc.issparse(nonlinear.jac_y):
                D = spc.lil_matrix((n_ineq, n_ineq))
                D.setdiag(sign)
                J = D*nonlinear.jac_y[ineq, :]
            else:
                J = np.multiply(nonlinear.jac_y[ineq, :],
                                sign[:, np.newaxis])
            return check_sparsity(J, sparse_jacobian)
        else:
            return np.empty((0, len(x)))

    if nonlinear.hess is None:
        hess = None
    else:
        def hess(x, v_eq=np.empty(0), v_ineq=np.empty(0)):
            hess = nonlinear.hess
            v = np.zeros(fun_len)
            if len(v_eq) > 0:
                v[eq] += v_eq
            if len(v_ineq) > 0:
                v[ineq[sign == 1]] += v_ineq[sign == 1]
                v[ineq[sign == -1]] -= v_ineq[sign == -1]
            return hess(x, v)

    if n_ineq == 0:
        enforce_feasibility = np.empty(0, dtype=bool)
    else:
        enforce_feasibility = nonlinear.enforce_feasibility[ineq]

    return CanonicalConstraint(n_ineq, constr_ineq, jac_ineq,
                               n_eq, constr_eq, jac_eq, hess,
                               enforce_feasibility)


def linear_to_canonical(linear, sparse_jacobian=None):
    return nonlinear_to_canonical(linear.to_nonlinear(sparse_jacobian))


def box_to_canonical(box, sparse_jacobian=None):
    return linear_to_canonical(box.to_linear(sparse_jacobian))


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


def empty_canonical_constraint():
    """Return empty CanonicalConstraint."""
    n_eq = 0
    n_ineq = 0

    def constr_eq(x):
        return np.empty(0)

    def constr_ineq(x):
        return np.empty(0)

    def jac_eq(x):
        return spc.csc_matrix(np.empty((0, len(x))))

    def jac_ineq(x):
        return spc.csc_matrix(np.empty((0, len(x))))

    enforce_feasibility = np.empty(0, dtype=bool)
    return CanonicalConstraint(n_ineq, constr_ineq, jac_ineq,
                               n_eq, constr_eq, jac_eq, None,
                               enforce_feasibility)


def generate_lagrangian_hessian(constraint, hess):
    """Lagrangian hessian"""

    # Concatenate Hessians
    def lagr_hess(x, v_eq=np.empty(0), v_ineq=np.empty(0)):
        n = len(x)
        hess_list = []
        if hess is not None:
            hess_list += [hess(x)]
        if constraint.hess is not None:
            hess_list += [constraint.hess(x, v_eq, v_ineq)]

        def matvec(p):
            result = np.zeros_like(p)
            for h in hess_list:
                result += h.dot(p)
            return result

        return spc.linalg.LinearOperator((n, n), matvec)

    return lagr_hess


def concatenate_canonical_constraints(constraints, sparse_jacobian=None):
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
        # Read sequentially all jacobians
        jac_eq_list = []
        sparse_jac_list = []
        for constr in constraints:
            J = constr.jac_eq(x)
            jac_eq_list += [J]
            sparse_jac_list += [spc.issparse(J)]
        # Use sparse matrix if any of the constraints
        # are sparse
        if sparse_jacobian is None:
            use_sparse = np.any(sparse_jac_list)
        else:
            use_sparse = sparse_jacobian
        # Convert all values to same format
        # this is done internally anyway and
        # it helps avoiding some odd behaviours
        for i in range(len(jac_eq_list)):
            if use_sparse:
                jac_eq_list[i] = spc.coo_matrix(jac_eq_list[i])
            elif sparse_jac_list[i]:
                jac_eq_list[i] = jac_eq_list[i].toarray()
        # Concatenate all
        if use_sparse:
            return spc.vstack(jac_eq_list, format="csc")
        else:
            return np.vstack(jac_eq_list)

    # Concatanate inequality constraints Jacobian matrices
    def jac_ineq(x):
        # Read sequentially all jacobians
        jac_ineq_list = []
        sparse_jac_list = []
        for constr in constraints:
            J = constr.jac_ineq(x)
            jac_ineq_list += [J]
            sparse_jac_list += [spc.issparse(J)]
        # Use sparse matrix if any of the constraints
        # are sparse in the case sparse_jacobian=None
        if sparse_jacobian is None:
            use_sparse = np.any(sparse_jac_list)
        else:
            use_sparse = sparse_jacobian
        # Convert all values to same format
        # this is done internally anyway and
        # it helps avoiding some odd behaviours
        for i in range(len(jac_ineq_list)):
            if use_sparse:
                jac_ineq_list[i] = spc.coo_matrix(jac_ineq_list[i])
            elif sparse_jac_list[i]:
                jac_ineq_list[i] = jac_ineq_list[i].toarray()
        # Concatenate all
        if use_sparse:
            return spc.vstack(jac_ineq_list, format="csc")
        else:
            return np.vstack(jac_ineq_list)

    # Concatenate Hessians
    def new_hess(x, v_eq=np.empty(0), v_ineq=np.empty(0)):
        n = len(x)
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

    # Concatenate feasible constraint list
    enforce_feasibility_list = [constr.enforce_feasibility
                                for constr in constraints]
    enforce_feasibility = np.hstack(enforce_feasibility_list)

    return CanonicalConstraint(n_ineq, constr_ineq, jac_ineq,
                               n_eq, constr_eq, jac_eq, new_hess,
                               enforce_feasibility)
