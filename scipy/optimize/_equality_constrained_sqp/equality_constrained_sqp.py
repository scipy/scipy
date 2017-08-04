"""Byrd-Omojokun Trust-Region SQP method."""

from __future__ import division, print_function, absolute_import
import scipy.sparse as spc
from .projections import projections
from .qp_subproblem import modified_dogleg, projected_cg, box_intersections
import numpy as np
from numpy.linalg import norm
from ..optimize import OptimizeResult

__all__ = ['equality_constrained_sqp']


def default_scaling(x):
    n, = np.shape(x)
    return spc.eye(n)


def equality_constrained_sqp(fun, grad, hess, constr, jac,
                             x0, stop_criteria,
                             trust_lb=None,
                             trust_ub=None,
                             initial_penalty=1.0,
                             initial_trust_radius=1.0,
                             scaling=default_scaling,
                             return_all=False):
    """Solve nonlinear equality-constrained problem using trust-region SQP.

    Solve problem:

        minimize fun(x)
        subject to: constr(x) = 0

    using Byrd-Omojokun Trust-Region SQP method.

    Parameters
    ----------
    fun : callable
        Objective function:
            fun(x) -> float
    grad : callable
        Gradient vector:
            grad(x) -> array_like, shape (n,)
    hess : callable
        Lagrangian hessian:
            hess(x, v) -> LinearOperator (or sparse matrix or ndarray), shape (n, n)
    constr : callable
        Equality constraint:
            constr(x) -> array_like, shape (m,)
    jac : callable
        Constraints Jacobian:
            jac(x) -> sparse matrix (or ndarray), shape (m, n)
    x0 : array_like, shape (n,)
        Starting point.
    v0 : array_like, shape (n,)
        Initial lagrange multipliers. By default uses least-squares lagrange
        multipliers.
    initial_trust_radius: float
        Initial trust-region radius. By defaut uses 1.
    trust_lb : array_like, shape (n,), optional
        Trust region lower bound.
    trust_ub : array_like, shape (n,), optional
        Trust region upper bound.
    stop_criteria: callable
        Functions that returns True when stop criteria is fulfilled:
            stop_criteria(state)
    initial_penalty : float
        Initial penalty for merit function.
    scaling : callable
        Function that return scaling used by the trust region:
            scaling(x) -> LinearOperator (or sparse matrix or ndarray), shape (n, n)
    return_all : bool, optional
        When ``true`` return the list of all vectors through the iterations.

    Returns
    -------
    x : array_like, shape (n,)
        Solution to the equality constrained problem.
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
            - allvecs : List containing all intermediary vectors (optional).
            - allmult : List containing all intermediary lagrange
                        multipliers (optional).

    Notes
    -----
    This algorithm is a variation of Byrd-Omojokun Trust-Region
    SQP method (described in [2]_ p.549 and in [3]_).Some details of
    this specific implementation are also inspired by [1]_.

    At each substep solve, using the projected CG method, the trust-region
    QP subproblem:

        minimize c.T d + 1/2 d.T H d
        subject to : A d + b = 0
                    ||d|| <= trust_radius
                    trust_lb <= d <= trust_ub

    and update the solution ``x += d``.

    References
    ----------
    .. [1] Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal.
           "An interior point algorithm for large-scale nonlinear
           programming." SIAM Journal on Optimization 9.4 (1999): 877-900.
    .. [2] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
           Second Edition (2006).
    .. [3] Lalee, Marucha, Jorge Nocedal, and Todd Plantenga. "On the
           implementation of an algorithm for large-scale equality
           constrained optimization." SIAM Journal on
           Optimization 8.3 (1998): 682-706.
    """
    PENALTY_FACTOR = 0.3  # Rho from formula (3.51), reference [1]_, p.891.
    LARGE_REDUCTION_RATIO = 0.9
    INTERMEDIARY_REDUCTION_RATIO = 0.3
    SUFFICIENT_REDUCTION_RATIO = 1e-8  # Eta from reference [1]_, p.892.
    TRUST_ENLARGEMENT_FACTOR_L = 7.0
    TRUST_ENLARGEMENT_FACTOR_S = 2.0
    MAX_TRUST_REDUCTION = 0.5
    MIN_TRUST_REDUCTION = 0.1
    SOC_THRESHOLD = 0.1
    TR_FACTOR = 0.8  # Zeta from formula (3.21), reference [1]_, p.885.
    BOX_FACTOR = 0.5

    n, = np.shape(x0)  # Number of parameters

    # Set default lower and upper bounds.
    if trust_lb is None:
        trust_lb = np.full(n, -np.inf)
    if trust_ub is None:
        trust_ub = np.full(n, np.inf)

    # Construct State structure
    state = OptimizeResult(niter=0, nfev=0, ngev=0,
                           ncev=0, njev=0, nhev=0)

    # Initial values
    x = np.copy(x0)
    trust_radius = initial_trust_radius
    penalty = initial_penalty
    # Compute Values
    f = fun(x)
    c = grad(x)
    b = constr(x)
    A = jac(x)
    S = scaling(x)
    # Get projections
    Z, LS, Y = projections(A)
    # Compute least-square lagrange multipliers
    v = -LS.dot(c)
    # Update state parameters
    state.optimality = norm(c + A.T.dot(v))
    state.constr_violation = norm(b)
    state.nfev += 1
    state.ngev += 1
    state.ncev += 1
    state.njev += 1
    state.x = x
    state.v = v
    state.fun = f
    state.grad = c
    state.constr = b
    state.jac = A
    state.trust_radius = trust_radius
    state.penalty = penalty
    # Store values
    if return_all:
        allvecs = [np.copy(x)]
        allmult = [np.copy(v)]

    compute_hess = True
    while not stop_criteria(state):
        # Compute Lagrangian Hessian
        if compute_hess:
            H = hess(x, v)
            state.nhev += 1

        # Normal Step - `dn`
        # minimize 1/2*||A dn + b||^2
        # subject to:
        # ||dn|| <= TR_FACTOR * trust_radius
        # BOX_FACTOR * lb <= dn <= BOX_FACTOR * ub.
        dn = modified_dogleg(A, Y, b,
                             TR_FACTOR*trust_radius,
                             BOX_FACTOR*trust_lb,
                             BOX_FACTOR*trust_ub)

        # Tangential Step - `dn`
        # Solve the QP problem:
        # minimize 1/2 dt.T H dt + dt.T (H dn + c)
        # subject to:
        # A dt = 0
        # ||dt|| <= sqrt(trust_radius**2 - ||dn||**2)
        # lb - dn <= dt <= ub - dn
        c_t = H.dot(dn) + c
        b_t = np.zeros_like(b)
        trust_radius_t = np.sqrt(trust_radius**2 - np.linalg.norm(dn)**2)
        lb_t = trust_lb - dn
        ub_t = trust_ub - dn
        dt, info_cg = projected_cg(H, c_t, Z, Y, b_t,
                                   trust_radius_t,
                                   lb_t, ub_t)

        # Compute update (normal + tangential steps).
        d = dn + dt

        # Compute second order model: 1/2 d H d + c.T d + f.
        quadratic_model = 1/2*(H.dot(d)).dot(d) + c.T.dot(d)
        # Compute linearized constraint: l = A d + b.
        linearized_constr = A.dot(d)+b
        # Compute new penalty parameter according to formula (3.52),
        # reference [1]_, p.891.
        vpred = norm(b) - norm(linearized_constr)
        # Guarantee `vpred` always positive,
        # regardless of roundoff errors.
        vpred = max(1e-16, vpred)
        previous_penalty = penalty
        if quadratic_model > 0:
            new_penalty = quadratic_model / ((1-PENALTY_FACTOR)*vpred)
            penalty = max(penalty, new_penalty)
        # Compute predicted reduction according to formula (3.52),
        # reference [1]_, p.891.
        predicted_reduction = -quadratic_model + penalty*vpred

        # Compute merit function at current point
        merit_function = f + penalty*norm(b)
        # Evaluate function and constraints at trial point
        x_next = x + S.dot(d)
        f_next = fun(x_next)
        b_next = constr(x_next)
        # Increment funcion evaluation counter
        state.nfev += 1
        state.ncev += 1
        # Compute merit function at trial point
        merit_function_next = f_next + penalty*norm(b_next)
        # Compute actual reduction according to formula (3.54),
        # reference [1]_, p.892.
        actual_reduction = merit_function - merit_function_next
        # Compute reduction ratio
        reduction_ratio = actual_reduction / predicted_reduction

        # Second order correction (SOC), reference [1]_, p.892.
        if reduction_ratio < SUFFICIENT_REDUCTION_RATIO and \
           norm(dn) <= SOC_THRESHOLD * norm(dt):
            # Compute second order correction
            y = -Y.dot(b_next)
            # Make sure increment is inside box constraints
            _, t, intersect = box_intersections(d, y, trust_lb, trust_ub)
            # Compute tentative point
            x_soc = x + S.dot(d + t*y)
            f_soc = fun(x_soc)
            b_soc = constr(x_soc)
            # Increment funcion evaluation counter
            state.nfev += 1
            state.ncev += 1
            # Recompute actual reduction
            merit_function_soc = f_soc + penalty*norm(b_soc)
            actual_reduction_soc = merit_function - merit_function_soc
            # Recompute reduction ratio
            reduction_ratio_soc = actual_reduction_soc / predicted_reduction
            if intersect and reduction_ratio_soc >= SUFFICIENT_REDUCTION_RATIO:
                x_next = x_soc
                f_next = f_soc
                b_next = b_soc
                reduction_ratio = reduction_ratio_soc

        # Readjust trust region step, formula (3.55), reference [1]_, p.892.
        if reduction_ratio >= LARGE_REDUCTION_RATIO:
            trust_radius = max(TRUST_ENLARGEMENT_FACTOR_L * norm(d),
                               trust_radius)
        elif reduction_ratio >= INTERMEDIARY_REDUCTION_RATIO:
            trust_radius = max(TRUST_ENLARGEMENT_FACTOR_S * norm(d),
                               trust_radius)
        # Reduce trust region step, according to reference [3]_, p.696.
        elif reduction_ratio < SUFFICIENT_REDUCTION_RATIO:
                trust_reduction = (1-SUFFICIENT_REDUCTION_RATIO) / (1-reduction_ratio)
                new_trust_radius = trust_reduction * norm(d)
                if new_trust_radius >= MAX_TRUST_REDUCTION * trust_radius:
                    trust_radius *= MAX_TRUST_REDUCTION
                elif new_trust_radius >= MIN_TRUST_REDUCTION * trust_radius:
                    trust_radius = new_trust_radius
                else:
                    trust_radius *= MIN_TRUST_REDUCTION

        # Update iteration
        state.niter += 1
        if reduction_ratio >= SUFFICIENT_REDUCTION_RATIO:
            x = x_next
            f = f_next
            c = grad(x)
            b = b_next
            A = jac(x)
            S = scaling(x)
            # Increment funcion evaluation counter
            state.ngev += 1
            state.njev += 1
            # Get projections
            Z, LS, Y = projections(A)
            # Compute least-square lagrange multipliers
            v = -LS.dot(c)
            # Set Flag
            compute_hess = True
            # Store state
            state.x = x
            state.v = v
            state.fun = f
            state.grad = c
            state.constr = b
            state.jac = A
            # Otimality values
            state.optimality = norm(c + A.T.dot(v))
            state.constr_violation = norm(b)
        else:
            penalty = previous_penalty
            compute_hess = False
        # Store values
        state.trust_radius = trust_radius
        state.penalty = penalty
        if return_all:
            allvecs.append(np.copy(x))
            allmult.append(np.copy(v))

    if return_all:
        state.allvecs = allvecs
        state.allmult = allmult

    return state
