from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc
from ._constraints import (NonlinearConstraint,
                           LinearConstraint,
                           BoxConstraint)

__all__ = ['CanonicalConstraint',
           'to_canonical',
           'lagrangian_hessian',
           'empty_canonical_constraint']


class CanonicalConstraint:
    """Canonical constraint

    Constraint of the form:
        c_ineq <= 0
        c_eq = 0
    for:
        c_ineq, c_eq = constr(x)
    """
    def __init__(self, n_vars, n_ineq, n_eq,
                 constr, jac, hess, sparse_jacobian,
                 enforce_feasibility,
                 x0, c_ineq0, c_eq0, J_ineq0, J_eq0):
        # Dimensions
        self.n_vars = n_vars
        self.n_ineq = n_ineq
        self.n_eq = n_eq
        # Objective function and constraints
        self.constr = constr
        self.jac = jac
        self.hess = hess
        # Use sparse jacobian flag
        self.sparse_jacobian = sparse_jacobian
        # Enforce feasibility for CanonicalConstraint. Should
        # be a list of booleans (and never a single boolean value,
        # as it is allowed for Box, Linear and Nonlinear constraints).
        self.enforce_feasibility = enforce_feasibility
        # Initial Values
        self.x0 = x0
        self.c_ineq0 = c_ineq0
        self.c_eq0 = c_eq0
        self.J_ineq0 = J_ineq0
        self.J_eq0 = J_eq0


def to_canonical(constraints):
    """Convert constraints or list of constraints to canonical format."""
    # Put ``constraints`` in list format whe needed
    if isinstance(constraints, (NonlinearConstraint,
                                LinearConstraint,
                                BoxConstraint,
                                CanonicalConstraint)):
        constraints = [constraints]
    if isinstance(constraints, (list, tuple, np.array)):
        # Converts all constraints to canonical format
        constraints_list = []
        for c in constraints:
            if isinstance(c, CanonicalConstraint):
                constraints_list += [c]
            elif isinstance(c, (NonlinearConstraint)):
                constraints_list += [_nonlinear_to_canonical(c)]
            elif isinstance(c, (LinearConstraint)):
                constraints_list += [_linear_to_canonical(c)]
            elif isinstance(c, (BoxConstraint)):
                constraints_list += [_box_to_canonical(c)]
            else:
                raise ValueError("Unknown Constraint type.")
        # Concatenate constraints
        if len(constraints_list) == 0:
            raise ValueError("Empty list.")
        elif len(constraints_list) == 1:
            constr = constraints_list[0]
        else:
            constr = _concatenate_canonical_constraints(constraints_list)
    else:
        raise ValueError("Unknown Constraint type.")

    return constr


def evaluated_to_canonical(constraints, list_c, list_J,
                           n_vars, n_eq, n_ineq, sparse_jacobian):
    n_constr = len(constraints)
    new_list_c = []
    new_list_J = []
    for i in range(n_constr):
        constr = constraints[i]
        c = list_c[i]
        J = list_J[i]
        eq, ineq, val_eq, val_ineq, sign, fun_len \
            = _parse_constraint(constr.kind)

        new_list_c += [_convert_constr(c, n_vars, n_eq, n_ineq,
                                       eq, ineq, val_eq, val_ineq,
                                       sign)]

        if constr.sparse_jacobian:
            new_list_J += [_convert_sparse_jac(J, n_vars, n_eq, n_ineq,
                                               eq, ineq, val_eq, val_ineq,
                                               sign)]
        else:
            new_list_J += [_convert_dense_jac(J, n_vars, n_eq, n_ineq,
                                              eq, ineq, val_eq, val_ineq,
                                              sign)]

    if sparse_jacobian:
        c_ineq, c_eq = _concatenate_constr(new_list_c)
        J_ineq, J_eq = _concatenate_sparse_jac(new_list_J)
    else:
        c_ineq, c_eq = _concatenate_constr(new_list_c)
        J_ineq, J_eq = _concatenate_dense_jac(new_list_J)

    return c_ineq, c_eq, J_ineq, J_eq


def lagrangian_hessian(constraint, hess):
    """Generate lagrangian hessian."""

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


def empty_canonical_constraint(x0, n_vars, sparse_jacobian=None):
    """Return empty CanonicalConstraint."""
    n_eq = 0
    n_ineq = 0
    empty_c = np.empty(0)
    if sparse_jacobian or (sparse_jacobian is None):
        empty_J = spc.csr_matrix(np.empty((0, n_vars)))
    else:
        empty_J = np.empty((0, n_vars))

    def constr(x):
        return empty_c, empty_c

    def jac(x):
        return empty_J, empty_J

    enforce_feasibility = np.empty(0, dtype=bool)
    return CanonicalConstraint(n_vars, n_ineq, n_eq,
                               constr, jac, None,
                               True, enforce_feasibility,
                               x0, empty_c, empty_c,
                               empty_J, empty_J)


# ************************************************************ #
# **********           Auxiliar Functions           ********** #
# ************************************************************ #
def _nonlinear_to_canonical(nonlinear):
    # Parse constraints
    eq, ineq, val_eq, val_ineq, sign, fun_len \
        = _parse_constraint(nonlinear.kind)
    # Get dimensions
    n_eq = len(eq)
    n_ineq = len(ineq)
    n_vars = nonlinear.n

    def new_constr(x):
        c = nonlinear.fun(x)
        return _convert_constr(c, n_vars, n_eq, n_ineq,
                               eq, ineq, val_eq, val_ineq,
                               sign)
    c_ineq0, c_eq0 = _convert_constr(nonlinear.f0, n_vars, n_eq, n_ineq,
                                     eq, ineq, val_eq, val_ineq,
                                     sign)

    if nonlinear.sparse_jacobian:
        def new_jac(x):
            J = nonlinear.jac(x)
            return _convert_sparse_jac(J, n_vars, n_eq, n_ineq,
                                       eq, ineq, val_eq, val_ineq,
                                       sign)
        J_ineq0, J_eq0 = _convert_sparse_jac(nonlinear.J0, n_vars, n_eq,
                                             n_ineq, eq, ineq, val_eq,
                                             val_ineq, sign)

    else:
        def new_jac(x):
            J = nonlinear.jac(x)
            return _convert_dense_jac(J, n_vars, n_eq, n_ineq,
                                      eq, ineq, val_eq, val_ineq,
                                      sign)
        J_ineq0, J_eq0 = _convert_dense_jac(nonlinear.J0, n_vars, n_eq,
                                            n_ineq, eq, ineq, val_eq,
                                            val_ineq, sign)

    if nonlinear.hess is None:
        new_hess = None
    else:
        def new_hess(x, v_eq=np.empty(0), v_ineq=np.empty(0)):
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

    return CanonicalConstraint(n_vars, n_ineq, n_eq,
                               new_constr, new_jac, new_hess,
                               nonlinear.sparse_jacobian,
                               enforce_feasibility, nonlinear.x0,
                               c_ineq0, c_eq0, J_ineq0, J_eq0)


def _linear_to_canonical(linear):
    return _nonlinear_to_canonical(linear.to_nonlinear())


def _box_to_canonical(box):
    return _linear_to_canonical(box.to_linear())


def _convert_constr(c, n_vars, n_eq, n_ineq,
                    eq, ineq, val_eq, val_ineq,
                    sign):
    # Empty constraint
    empty = np.empty((0,))
    # Return equality and inequalit constraints
    c_eq = c[eq] - val_eq if n_eq > 0 else empty
    c_ineq = sign*(c[ineq] - val_ineq) if n_ineq > 0 else empty
    return c_ineq, c_eq


def _convert_sparse_jac(J, n_vars, n_eq, n_ineq,
                        eq, ineq, val_eq, val_ineq,
                        sign):
    # Empty jacobian
    empty = spc.csr_matrix(np.empty((0, n_vars)))
    # Compute equality and inequality Jacobian matrices
    J_eq = J[eq, :] if n_eq > 0 else empty
    if n_ineq > 0:
        D = spc.lil_matrix((n_ineq, n_ineq))
        D.setdiag(sign)
        J_ineq = D*J[ineq, :]
    else:
        J_ineq = empty
    # Return Jacobian matrices
    return J_ineq, J_eq


def _convert_dense_jac(J, n_vars, n_eq, n_ineq,
                       eq, ineq, val_eq, val_ineq,
                       sign):
    # Empty jacobian
    empty = np.empty((0, n_vars))
    # Compute equality and inequality Jacobian matrices
    J_eq = J[eq, :] if n_eq > 0 else empty
    if n_ineq > 0:
        J_ineq = np.multiply(J[ineq, :], sign[:, np.newaxis])
    else:
        J_ineq = empty
    # Return Jacobian matrices
    return J_ineq, J_eq


def _parse_constraint(kind):
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
        ineq = np.empty(0, dtype=int)
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
        raise RuntimeError("Never be here.")

    return eq, ineq, val_eq, val_ineq, sign, fun_len


def _concatenate_canonical_constraints(constraints,
                                       sparse_jacobian=None):
    """Concatenate sequence of CanonicalConstraint's."""
    # Compute number of constraints
    n_eq = 0
    n_ineq = 0
    for constr in constraints:
        n_eq += constr.n_eq
        n_ineq += constr.n_ineq

    # Get n_vars
    n_vars = 0
    x0 = None
    for constr in constraints:
        if n_vars == 0:
            n_vars = constr.n_vars
            x0 = constr.x0
        if n_vars != constr.n_vars:
            raise RuntimeError("Unmatching constraint number of arguments.")
        if not np.array_equal(x0, constr.x0):
            raise RuntimeError("Unmatching initial point.")

    # Concatenate constraints
    def new_constr(x):
        constr_list = [constr.constr(x) for constr in constraints]
        return _concatenate_constr(constr_list)
    constr0_list = [(constr.c_ineq0, constr.c_eq0) for constr in constraints]
    c_ineq0, c_eq0 = _concatenate_constr(constr0_list)

    # Use sparse if any of the matrices are sparse
    use_sparse = np.any([constr.sparse_jacobian for constr in constraints])

    if use_sparse:
        def new_jac(x):
            jac_list = [constr.jac(x) for constr in constraints]
            return _concatenate_sparse_jac(jac_list)
        jac0_list = [(constr.J_ineq0, constr.J_eq0) for constr in constraints]
        J_ineq0, J_eq0 = _concatenate_sparse_jac(jac0_list)

    else:
        def new_jac(x):
            jac_list = [constr.jac(x) for constr in constraints]
            return _concatenate_dense_jac(jac_list)
        jac0_list = [(constr.J_ineq0, constr.J_eq0) for constr in constraints]
        J_ineq0, J_eq0 = _concatenate_dense_jac(jac0_list)

    # Concatenate Hessians
    def new_hess(x, v_eq=np.empty(0), v_ineq=np.empty(0)):
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

        return spc.linalg.LinearOperator((n_vars, n_vars), matvec)

    # Concatenate feasible constraint list
    enforce_feasibility_list = [constr.enforce_feasibility
                                for constr in constraints]
    enforce_feasibility = np.hstack(enforce_feasibility_list)

    return CanonicalConstraint(n_vars, n_ineq, n_eq, new_constr,
                               new_jac, new_hess, use_sparse,
                               enforce_feasibility, x0, c_ineq0,
                               c_eq0, J_ineq0, J_eq0)


def _concatenate_constr(constr_list):
    c_ineq = np.hstack([constr[0] for constr in constr_list])
    c_eq = np.hstack([constr[1] for constr in constr_list])
    return c_ineq, c_eq


def _concatenate_sparse_jac(jac_list):
    jac_ineq_list = []
    jac_eq_list = []
    for jac_tuple in jac_list:
        J_ineq, J_eq = jac_tuple
        jac_ineq_list += [spc.csr_matrix(J_ineq)]
        jac_eq_list += [spc.csr_matrix(J_eq)]
    # Concatenate all
    J_ineq = spc.vstack(jac_ineq_list, format="csr")
    J_eq = spc.vstack(jac_eq_list, format="csr")
    # Return
    return J_ineq, J_eq


def _concatenate_dense_jac(jac_list):
    # Read sequentially all jacobians.
    # Convert all values to numpy arrays.
    jac_ineq_list = []
    jac_eq_list = []
    for jac_tuple in jac_list:
        J_ineq, J_eq = jac_tuple
        if spc.issparse(J_ineq):
            jac_ineq_list += [J_ineq.toarray()]
        else:
            jac_ineq_list += [np.atleast_2d(J_ineq)]
        if spc.issparse(J_eq):
            jac_eq_list += [J_eq.toarray()]
        else:
            jac_eq_list += [np.atleast_2d(J_eq)]
    # Concatenate all
    J_ineq = np.vstack(jac_ineq_list)
    J_eq = np.vstack(jac_eq_list)
    # Return
    return J_ineq, J_eq
