"""Trust Region Reflective algorithm for least-squares optimization."""

from __future__ import division

import numpy as np
from numpy.linalg import norm
from scipy.linalg import svd

from .optimize import OptimizeResult
from ._lsq_bounds  import (step_size_to_bound, make_strictly_feasible,
                           find_active_constraints, scaling_vector)
from ._lsq_trust_region import intersect_trust_region, solve_lsq_trust_region


def minimize_quadratic(a, b, l, u):
    """Minimize a 1-d quadratic function subject to bounds.

    The free term is omitted, that is we consider y = a * t**2 + b * t.

    Returns
    -------
    t : float
        The minimum point.
    y : float
        The minimum value.
    """
    t = np.array([l, u])
    if a != 0:
        extremum = -0.5 * b / a
        if l <= extremum <= u:
            t = np.hstack((t, extremum))
    y = a * t**2 + b * t
    i = np.argmin(y)
    return t[i], y[i]


def build_1d_quadratic_function(J, diag, g, s, s0=None):
    """Compute coefficients of a 1-d quadratic function for the line search
    from the multidimensional quadratic function.

    The function is given as follows:
    ``f(t) = 0.5 * (s0 + t*s).T * (J.T*J + diag) * (s0 + t*s) +
             g.T * (s0 + t*s)``.

    Parameters
    ----------
    J : ndarray, shape (m, n)
        Jacobian matrix.
    diag : ndarray, shape (n,)
        Addition diagonal term in a quadratic part.
    g : ndarray, shape (n,)
        Gradient, defines a linear term.
    s : ndarray, shape (n,)
        Direction of search.
    s0 : None or ndarray with shape (n,), optional
        Initial point. If None assumed to be 0.

    Returns
    ------
    a : float
        Coefficient for t**2.
    b : float
        Coefficient for t.

    Notes
    -----
    The free term "c" is not returned as it is not usually required.
    """
    v = J.dot(s)
    a = 0.5 * (np.dot(v, v) + np.dot(s * diag, s))
    b = np.dot(g, s)
    if s0 is not None:
        u = J.dot(s0)
        b += np.dot(u, v) + np.dot(s0 * diag, s)

    return a, b


def evaluate_quadratic_function(J, diag, g, steps):
    """Compute values of a quadratic function arising in least-squares.

    The function is 0.5 * s.T * (J.T * J + diag) * s + g.T * s.

    Parameters
    ----------
    J : ndarray, shape (m, n)
        Jacobian matrix.
    diag : ndarray, shape (n,)
        Additional diagonal term.
    f : ndarray, shape (m,)
        Vector of residuals.
    steps : ndarray, shape (k, m)
        Array containing k steps as rows.

    Returns
    -------
    values : ndarray, shape (k,)
        Array containing k values of the function.
    """
    Js = J.dot(steps.T)
    return 0.5 * (np.sum(Js**2, axis=0) +
                  np.sum(diag * steps**2, axis=1)) + np.dot(steps, g)


def find_reflected_step(x, J_h, diag_h, g_h, p, p_h, d, Delta, l, u, theta):
    """Find a single reflection step for Trust Region Reflective algorithm.

    Also corrects the initial step p/p_h. This function must be called only
    if x + p is not within the bounds.
    """
    # Use term "stride" for scalar step length.
    p_stride, hits = step_size_to_bound(x, p, l, u)

    # Compute the reflected direction.
    r_h = np.copy(p_h)
    r_h[hits.astype(bool)] *= -1
    r = d * r_h

    # Restrict the p step to exactly the bound.
    p *= p_stride
    p_h *= p_stride
    x_on_bound = x + p

    # Reflected direction will cross first either feasible region or trust
    # region boundary.
    _, to_tr = intersect_trust_region(p_h, r_h, Delta)
    to_bound, _ = step_size_to_bound(x_on_bound, r, l, u)
    to_bound *= theta  # Stay strictly interior.

    r_stride_u = min(to_bound, to_tr)

    # We want a reflected step be at the same theta distance from the bound,
    # so we introduce a lower bound on the allowed stride.
    # The formula below is correct as p_h and r_h has the same norm.

    if r_stride_u > 0:
        r_stride_l = (1 - theta) * p_stride / r_stride_u
    else:
        r_stride_l = -1

    # Check if no reflection step is available.
    if r_stride_l <= r_stride_u:
        a, b = build_1d_quadratic_function(J_h, diag_h, g_h, r_h, s0=p_h)
        r_stride, _ = minimize_quadratic(a, b, r_stride_l, r_stride_u)
        r_h = p_h + r_h * r_stride
    else:
        r_h = None

    # Now we want to correct p_h to make it strictly interior.
    p_h *= theta

    # If no reflection step just return p_h for convenience.
    if r_h is None:
        return p_h, p_h
    else:
        return p_h, r_h


def find_gradient_step(x, J_h, diag_h, g_h, d, Delta, l, u, theta):
    """Find a minimizer of a quadratic model along the scaled gradient."""
    to_bound, _ = step_size_to_bound(x, -g_h * d, l, u)
    to_bound *= theta

    to_tr = Delta / norm(g_h)
    g_stride = min(to_bound, to_tr)

    a, b = build_1d_quadratic_function(J_h, diag_h, g_h, -g_h)
    g_stride, _ = minimize_quadratic(a, b, 0.0, g_stride)

    return -g_stride * g_h


def trf(fun, jac, x0, lb, ub, ftol, xtol, gtol, max_nfev, scaling):
    """Minimize the sum of squares of nonlinear functions subject to bounds on
    independent variables by Trust Region Reflective algorithm.

    Options
    -------
    ftol : float
        The optimization process is stopped when ``dF < ftol * F`` and
        dF_actual / dF_predicted > 0.25, where F is the objective function
        value (the sum of squares), dF_actual is its change in the last
        iteration, dF_predicted is predicted change from a local quadratic
        model.
    xtol : float, optional
        The optimization process is stopped when
        ``norm(dx) < xtol * max(EPS**0.5, norm(x))``, where dx is a step taken
        in the last iteration and EPS is machine epsilon.
    gtol : float, optional
        The optimization process is stopped when
        ``norm(g_scaled, ord=np.inf) < gtol``, where g_scaled is properly
        scaled gradient to account for the presence of bounds.
        The scaling imposed by `scaling` parameter is not considered.
    max_nfev : None or int, optional
        Maximum number of function evaluations before the termination.
        If None (default), it is assigned to 100 * n.
    """
    EPS = np.finfo(float).eps

    # We need strictly feasible guess to start with
    x = make_strictly_feasible(x0, lb, ub, rstep=1e-10)

    f = fun(x)
    nfev = 1

    J = jac(x, f)
    njev = 1

    if f.shape[0] != J.shape[0]:
        raise RuntimeError("Inconsistent dimensions between the returns of "
                           "`fun` and `jac` on the first iteration.")

    g = J.T.dot(f)
    m, n = J.shape

    if scaling == 'jac':
        J_norm = np.sum(J**2, axis=0)**0.5
        J_norm[J_norm == 0] = 1
        scale = 1 / J_norm
    else:
        scale = 1 / scaling

    v, jv = scaling_vector(x, g, lb, ub)
    Delta = norm(x0 / (scale * v**0.5))
    if Delta == 0:
        Delta = 1.0

    J_augmented = np.empty((m + n, n))
    f_augmented = np.zeros((m + n))

    obj_value = np.dot(f, f)
    alpha = 0.0  # "Levenberg-Marquardt" parameter

    if max_nfev is None:
        max_nfev = x0.size * 100

    termination_status = None
    while nfev < max_nfev:
        if scaling == 'jac':
            J_norm = np.sum(J**2, axis=0)**0.5
            with np.errstate(divide='ignore'):
                scale = np.minimum(scale, 1 / J_norm)

        g = J.T.dot(f)

        # Compute Coleman-Li scaling parameters and "hat" variables.
        v, jv = scaling_vector(x, g, lb, ub)
        d = v**0.5 * scale
        g_h = d * g
        diag_h = g * jv * scale**2

        g_norm = norm(g * v, ord=np.inf)
        if g_norm < gtol:
            termination_status = 1

        if termination_status is not None:
            active_mask = find_active_constraints(x, lb, ub, rtol=xtol)
            return OptimizeResult(
                x=x, fun=f, jac=J, obj_value=obj_value, optimality=g_norm,
                active_mask=active_mask, nfev=nfev, njev=njev,
                status=termination_status, x_covariance=None)

        # Jacobian in "hat" space.
        J_h = J * d

        # J_augmented is used to solve a trust-region subproblem with
        # diagonal term diag_h.
        J_augmented[:m] = J_h
        J_augmented[m:] = np.diag(diag_h**0.5)
        f_augmented[:m] = f

        U, s, V = svd(J_augmented, full_matrices=False)
        V = V.T
        uf = U.T.dot(f_augmented)

        # theta controls step back size from the bounds.
        theta = max(0.995, 1 - g_norm)
        actual_reduction = -1

        # In the following: p - trust-region solution, r - reflected solution,
        # c - minimizer along the scaled gradient, _h means the variable
        # is computed in "hat" space.
        while actual_reduction <= 0 and nfev < max_nfev:
            p_h, alpha, n_iter = solve_lsq_trust_region(
                n, m, uf, s, V, Delta, initial_alpha=alpha)
            p = d * p_h

            to_bound, _ = step_size_to_bound(x, p, lb, ub)
            if to_bound >= 1:  # Trust region step is feasible.
                # Still step back from the bounds
                p_h *= min(theta * to_bound, 1)
                steps_h = np.atleast_2d(p_h)
            else:  # Otherwise consider a reflected and gradient steps.
                p_h, r_h = find_reflected_step(x, J_h, diag_h, g_h, p, p_h,
                                               d, Delta, lb, ub, theta)
                c_h = find_gradient_step(x, J_h, diag_h, g_h,
                                         d, Delta, lb, ub, theta)
                steps_h = np.array([p_h, r_h, c_h])

            qp_values = evaluate_quadratic_function(J_h, diag_h, g_h, steps_h)
            min_index = np.argmin(qp_values)
            step_h = steps_h[min_index]

            # qp-values are negative, also need to double it
            predicted_reduction = -2 * qp_values[min_index]

            step = d * step_h
            x_new = make_strictly_feasible(x + step, lb, ub)

            f_new = fun(x_new)
            nfev += 1

            # Usual trust-region step quality estimation.
            obj_value_new = np.dot(f_new, f_new)
            actual_reduction = obj_value - obj_value_new
            # Correction term is specific to the algorithm,
            # vanishes in unbounded case.
            correction = np.dot(step_h * diag_h, step_h)

            if predicted_reduction > 0:
                ratio = (actual_reduction - correction) / predicted_reduction
            else:
                ratio = 0

            if ratio < 0.25:
                Delta_new = 0.25 * norm(step_h)
                alpha *= Delta / Delta_new
                Delta = Delta_new
            elif ratio > 0.75 and norm(step_h) > 0.95 * Delta:
                Delta *= 2.0
                alpha *= 0.5

            ftol_satisfied = (abs(actual_reduction) < ftol * obj_value and
                              ratio > 0.25)
            xtol_satisfied = norm(step) < xtol * max(EPS**0.5, norm(x))

            if ftol_satisfied and xtol_satisfied:
                termination_status = 4
            elif ftol_satisfied:
                termination_status = 2
            elif xtol_satisfied:
                termination_status = 3
            if termination_status is not None:
                break

        if actual_reduction > 0:
            x = x_new
            f = f_new
            obj_value = obj_value_new

            J = jac(x, f)
            njev += 1

    active_mask = find_active_constraints(x, lb, ub, rtol=xtol)
    return OptimizeResult(
        x=x, fun=f, jac=J, obj_value=obj_value, optimality=g_norm,
        active_mask=active_mask, nfev=nfev, njev=njev, status=0,
        x_covariance=None)
