from __future__ import division

from warnings import warn

import numpy as np
from numpy.linalg import lstsq, norm

from . import OptimizeResult
from ..sparse import issparse
from ..sparse.linalg import LinearOperator, aslinearoperator, lsmr
from ._lsq_bounds import step_size_to_bound, in_bounds


def lsmr_linear_operator(Jop, active_set):
    m, n = Jop.shape

    def matvec(x):
        x_free = x.copy()
        x_free[active_set] = 0
        return Jop.matvec(x)

    def rmatvec(x):
        r = Jop.rmatvec(x)
        r[active_set] = 0
        return r

    return LinearOperator((m, n), matvec=matvec, rmatvec=rmatvec, dtype=float)


def find_intersection(x, tr_bounds, lb, ub):
    """Find intersection of trust-region bounds and initial bounds.

    Returns
    -------
    lb_total, ub_total : ndarray with shape of x
        Lower and upper bounds of the intersection region.
    orig_l, orig_u : ndarray of bool with shape of x
        True means that an original bound is taken as a corresponding bound
        in the intersection region.
    tr_l, tr_u : ndarray of bool with shape of x
        True means that a trust-region bound is taken as a corresponding bound
        in the intersection region.
    """
    lb_centered = lb - x
    ub_centered = ub - x

    lb_total = np.maximum(lb_centered, -tr_bounds)
    ub_total = np.minimum(ub_centered, tr_bounds)

    orig_l = np.equal(lb_total, lb_centered)
    orig_u = np.equal(ub_total, ub_centered)

    tr_l = np.equal(lb_total, -tr_bounds)
    tr_u = np.equal(ub_total, tr_bounds)

    return lb_total, ub_total, orig_l, orig_u, tr_l, tr_u


def dogleg_step(x, cauchy_step, newton_step, tr_bounds, lb, ub):
    """Find dogleg step in rectangular region.

    Returns
    -------
    step : ndarray, shape (n,)
        Computed dogleg step.
    bound_hits : ndarray of int, shape (n,)
        Each component shows whether a corresponding variable hits the
        initial bound after the step is taken:

            *  0 - a variable doesn't hit the bound.
            * -1 - lower bound is hit.
            *  1 - upper bound is hit.
    tr_hit : bool
        Whether the step hit the boundary of the trust-region.
    """
    lb_total, ub_total, orig_l, orig_u, tr_l, tr_u = find_intersection(
        x, tr_bounds, lb, ub
    )
    bound_hits = np.zeros_like(x, dtype=int)

    if in_bounds(newton_step, lb_total, ub_total):
        return newton_step, bound_hits, False

    if not in_bounds(cauchy_step, lb_total, ub_total):
        step_size, _ = step_size_to_bound(
            np.zeros_like(cauchy_step), cauchy_step, lb_total, ub_total)
        cauchy_step = step_size * cauchy_step  # Don't want to modify inplace.
        # The classical dogleg algorithm would stop here, but in a rectangular
        # region it makes sense to try to improve constrained cauchy step.
        # Thus the code after this "if" is always executed.

    step_diff = newton_step - cauchy_step
    step_size, hits = step_size_to_bound(cauchy_step, step_diff,
                                         lb_total, ub_total)
    bound_hits[(hits < 0) & orig_l] = -1
    bound_hits[(hits > 0) & orig_u] = 1
    tr_hit = np.any((hits < 0) & tr_l | (hits > 0) & tr_u)

    return cauchy_step + step_size * step_diff, bound_hits, tr_hit


def constrained_cauchy_step(x, cauchy_step, tr_bounds, l, u):
    """Find constrained Cauchy step.

    Returns are the same as in dogleg_step function.
    """
    lb_total, ub_total, orig_l, orig_u, tr_l, tr_u = find_intersection(
        x, tr_bounds, l, u
    )
    bound_hits = np.zeros_like(x, dtype=int)

    if in_bounds(cauchy_step, lb_total, ub_total):
        return cauchy_step, bound_hits, False

    step_size, hits = step_size_to_bound(
        np.zeros_like(cauchy_step), cauchy_step, lb_total, ub_total)

    bound_hits[(hits < 0) & orig_l] = -1
    bound_hits[(hits > 0) & orig_u] = 1
    tr_hit = np.any((hits < 0) & tr_l | (hits > 0) & tr_u)

    return step_size * cauchy_step, bound_hits, tr_hit


def dogbox(fun, jac, x0, lb, ub, ftol, xtol, gtol, max_nfev, scaling,
           tr_solver, tr_options):
    """Minimize the sum of squares of nonlinear functions subject to bounds on
    independent variables by dogleg method applied to a rectangular trust
    region.

    Options
    -------
    ftol : float
        The optimization process is stopped when ``dF < ftol * F`` and
        dF_actual / dF_predicted > 0.25, where F is the objective function
        value (the sum of squares), dF_actual is its change in the last
        iteration, dF_predicted is predicted change from a local quadratic
        model.
    xtol : float
        The optimization process is stopped when
        ``Delta < xtol * max(EPS**0.5, norm(scaled_x))``, where Delta is a
        trust-region radius, scaled_x is a scaled value of x according
        to `scaling` parameter, EPS is machine epsilon.
    gtol : float
        The optimization process is stopped when
        ``norm(g_free, ord=np.inf) < gtol``, where g_free is the gradient
        with respect to the variables which aren't in the optimal state
        on the boundary. If all variables reach optimum on the boundary,
        then g_free is effectively assigned to zero and the algorithm
        terminates.
    max_nfev : None or int
        Maximum number of function evaluations before the termination.
        If None (default), it is assigned to 100 * n.
    """
    EPS = np.finfo(float).eps

    f = fun(x0)
    nfev = 1

    J = jac(x0, f)
    njev = 1
    if tr_solver is None:
        if issparse(J):
            tr_solver = 'lsmr'
        else:
            tr_solver = 'exact'
    elif tr_solver == 'exact' and issparse(J):
        warn("Sparse Jacobian will be converted to dense for tr_solver=exact, "
             "consider using 'lsmr' solver or return dense Jacobian.")
        J = J.toarray()

    if scaling == 'jac':
        scale = np.sum(J**2, axis=0)**0.5
        scale[scale == 0] = 1
    else:
        scale = scaling

    Delta = norm(x0 * scale, ord=np.inf)
    if Delta == 0:
        Delta = 1.0

    on_bound = np.zeros_like(x0, dtype=int)
    on_bound[np.equal(x0, lb)] = -1
    on_bound[np.equal(x0, ub)] = 1

    x = x0.copy()
    step = np.empty_like(x0)
    obj_value = np.dot(f, f)

    if max_nfev is None:
        max_nfev = x0.size * 100

    termination_status = None
    while nfev < max_nfev:
        if scaling == 'jac':
            scale = np.maximum(scale, np.sum(J**2, axis=0)**0.5)

        g = J.T.dot(f)

        active_set = on_bound * g < 0
        free_set = ~active_set

        g_free = g[free_set]
        if np.all(active_set):
            g_norm = 0.0
            termination_status = 1
        else:
            g_norm = norm(g_free, ord=np.inf)
            if g_norm < gtol:
                termination_status = 1

        if termination_status is not None:
            return OptimizeResult(
                x=x, fun=f, jac=J, obj_value=obj_value, optimality=g_norm,
                active_mask=on_bound, nfev=nfev, njev=njev,
                status=termination_status, x_covariance=None)

        x_free = x[free_set]
        l_free = lb[free_set]
        u_free = ub[free_set]
        scale_free = scale[free_set]

        # Compute (Gauss-)Newton and Cauchy steps
        if tr_solver == 'exact':
            J_free = J[:, free_set]
            newton_step = lstsq(J_free, -f)[0]
            Jg = J_free.dot(g_free)
        elif tr_solver == 'lsmr':
            Jop = aslinearoperator(J)
            lsmr_op = lsmr_linear_operator(Jop, active_set)
            newton_step = -lsmr(lsmr_op, f, **tr_options)[0][free_set]
            g[active_set] = 0
            Jg = Jop.matvec(g)
        cauchy_step = -np.dot(g_free, g_free) / np.dot(Jg, Jg) * g_free

        actual_reduction = -1.0
        while actual_reduction <= 0 and nfev < max_nfev:
            tr_bounds = Delta / scale_free

            step_free, on_bound_free, tr_hit = dogleg_step(
                x_free, cauchy_step, newton_step, tr_bounds, l_free, u_free)

            step.fill(0.0)
            step[free_set] = step_free

            if tr_solver == 'exact':
                Js = J_free.dot(step_free)
            elif tr_solver == 'lsmr':
                Js = Jop.matvec(step)

            predicted_reduction = -np.dot(Js, Js) - 2 * np.dot(Js, f)

            # In (nearly) rank deficient case Newton (and thus dogleg) step
            # can be inadequate, in this case use (constrained) Cauchy step.
            if predicted_reduction <= 0:
                step_free, on_bound_free, tr_hit = constrained_cauchy_step(
                    x_free, cauchy_step, tr_bounds, l_free, u_free)

                step.fill(0.0)
                step[free_set] = step_free

                if tr_solver == 'exact':
                    Js = J_free.dot(step_free)
                elif tr_solver == 'lsmr':
                    Js = Jop.matvec(step)

                predicted_reduction = -np.dot(Js, Js) - 2 * np.dot(Js, f)

            x_new = x + step
            f_new = fun(x_new)
            nfev += 1

            # Usual trust-region step quality estimation.
            obj_value_new = np.dot(f_new, f_new)
            actual_reduction = obj_value - obj_value_new

            if predicted_reduction > 0:
                ratio = actual_reduction / predicted_reduction
            else:
                ratio = 0

            if ratio < 0.25:
                Delta = 0.25 * norm(step * scale, ord=np.inf)
            elif ratio > 0.75 and tr_hit:
                Delta *= 2.0

            ftol_satisfied = (abs(actual_reduction) < ftol * obj_value and
                              ratio > 0.25)
            xtol_satisfied = Delta < xtol * max(EPS**0.5,
                                                norm(x * scale, ord=np.inf))

            if ftol_satisfied and xtol_satisfied:
                termination_status = 4
            elif ftol_satisfied:
                termination_status = 2
            elif xtol_satisfied:
                termination_status = 3
            if termination_status is not None:
                break

        if actual_reduction > 0:
            on_bound[free_set] = on_bound_free

            x = x_new
            # Set variables exactly at the boundary.
            mask = on_bound == -1
            x[mask] = lb[mask]
            mask = on_bound == 1
            x[mask] = ub[mask]

            f = f_new
            obj_value = obj_value_new

            J = jac(x, f)
            njev += 1
            if tr_solver == 'exact' and issparse(J):
                J = J.toarray()

    return OptimizeResult(
        x=x, fun=f, jac=J, obj_value=obj_value, optimality=g_norm,
        active_mask=on_bound, nfev=nfev, njev=njev, status=0,
        x_covariance=None)
