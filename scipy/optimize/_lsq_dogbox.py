"""
dogleg algorithm with rectangular trust regions for least-squares minimization.
"""

from __future__ import division

import numpy as np
from numpy.linalg import lstsq, norm

from . import OptimizeResult
from ..sparse import issparse
from ..sparse.linalg import LinearOperator, aslinearoperator, lsmr
from ._lsq_common import (step_size_to_bound, in_bounds,
                          print_header, print_iteration)


def lsmr_linear_operator(Jop, d, active_set):
    m, n = Jop.shape

    def matvec(x):
        x_free = x.copy()
        x_free[active_set] = 0
        return Jop.matvec(x * d)

    def rmatvec(x):
        r = d * Jop.rmatvec(x)
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


def dogbox(fun, jac, x0, f0, J0, lb, ub, ftol, xtol, gtol, max_nfev, scaling,
           tr_solver, tr_options, verbose):
    f = f0
    nfev = 1

    J = J0
    njev = 1

    if scaling == 'jac':
        if issparse(J):
            scale = np.asarray(J.power(2).sum(axis=0)).ravel()**0.5
        else:
            scale = np.sum(J**2, axis=0)**0.5
        scale[scale == 0] = 1
    else:
        scale = scaling
    scale_inv = 1 / scale

    Delta = norm(x0 * scale, ord=np.inf)
    if Delta == 0:
        Delta = 1.0

    on_bound = np.zeros_like(x0, dtype=int)
    on_bound[np.equal(x0, lb)] = -1
    on_bound[np.equal(x0, ub)] = 1

    x = x0.copy()
    step = np.empty_like(x0)
    obj_value = 0.5 * np.dot(f, f)

    if max_nfev is None:
        max_nfev = x0.size * 100

    termination_status = None
    iteration = 0
    step_norm = None
    actual_reduction = None

    if verbose == 2:
        print_header()

    while nfev < max_nfev:
        if scaling == 'jac':
            if issparse(J):
                new_scale = np.asarray(J.power(2).sum(axis=0)).ravel()**0.5
            else:
                new_scale = np.sum(J**2, axis=0)**0.5
            scale = np.maximum(scale, new_scale)
            scale_inv = 1 / scale

        if isinstance(J, LinearOperator):
            g = J.rmatvec(f)
        else:
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

        if verbose == 2:
            print_iteration(iteration, nfev, obj_value, g_norm,
                            step_norm, actual_reduction)

        if termination_status is not None:
            return OptimizeResult(
                x=x, fun=f, jac=J, obj_value=obj_value, optimality=g_norm,
                active_mask=on_bound, nfev=nfev, njev=njev,
                status=termination_status)

        x_free = x[free_set]
        l_free = lb[free_set]
        u_free = ub[free_set]
        scale_inv_free = scale_inv[free_set]

        # Compute (Gauss-)Newton and Cauchy steps
        if tr_solver == 'exact':
            J_free = J[:, free_set]
            newton_step = lstsq(J_free, -f)[0]
            Jg = J_free.dot(g_free)
        elif tr_solver == 'lsmr':
            Jop = aslinearoperator(J)
            # Here we compute lsmr step in scaled variables and then
            # transform back to normal variables, if lsmr would give exact lsq
            # solution this would be equivalent to not doing any
            # transformations, but from experience it's better this way.
            lsmr_op = lsmr_linear_operator(Jop, scale_inv, active_set)
            newton_step = -lsmr(lsmr_op, f, **tr_options)[0][free_set]
            newton_step *= scale_inv_free
            g[active_set] = 0
            Jg = Jop.matvec(g)
        cauchy_step = -np.dot(g_free, g_free) / np.dot(Jg, Jg) * g_free

        actual_reduction = -1.0
        while actual_reduction <= 0 and nfev < max_nfev:
            tr_bounds = Delta * scale_inv_free

            step_free, on_bound_free, tr_hit = dogleg_step(
                x_free, cauchy_step, newton_step, tr_bounds, l_free, u_free)

            step.fill(0.0)
            step[free_set] = step_free

            if tr_solver == 'exact':
                Js = J_free.dot(step_free)
            elif tr_solver == 'lsmr':
                Js = Jop.matvec(step)

            predicted_reduction = -0.5 * np.dot(Js, Js) - np.dot(Js, f)

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

                predicted_reduction = -0.5 * np.dot(Js, Js) - np.dot(Js, f)

            x_new = x + step
            f_new = fun(x_new)
            nfev += 1

            # Usual trust-region step quality estimation.
            obj_value_new = 0.5 * np.dot(f_new, f_new)
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

            step_norm = norm(step)
            xtol_satisfied = step_norm < xtol * (xtol + norm(x))

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
        iteration += 1

    return OptimizeResult(
        x=x, fun=f, jac=J, obj_value=obj_value, optimality=g_norm,
        active_mask=on_bound, nfev=nfev, njev=njev, status=0)
