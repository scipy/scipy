"""Fifth order implicit Runge-Kutta method of Radau family."""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from .common import (EPS, norm, get_active_events, handle_events,
                     select_initial_step)

S6 = 6 ** 0.5

# Butcher tableau. A is not used directly, see below.
C = np.array([(4 - S6) / 10, (4 + S6) / 10, 1])
E = np.array([-13 - 7 * S6, -13 + 7 * S6, -1]) / 3

# Eigendecomposition of A is done: A = T L T**-1. There is 1 real eigenvalue
# and a complex conjugate pair. They are written below.
MU_REAL = 3 + 3 ** (2 / 3) - 3 ** (1 / 3)
MU_COMPLEX = 3 + 0.5 * (3 ** (1 / 3) - 3 ** (2 / 3)) - 0.5j * (
    3 ** (5 / 6) + 3 ** (7 / 6))

# These are transformation matrices.
T = np.array([
    [0.09443876248897524, -0.14125529502095421, 0.03002919410514742],
    [0.25021312296533332, 0.20412935229379994, -0.38294211275726192],
    [1, 1, 0]])
TI = np.array([
    [4.17871859155190428, 0.32768282076106237, 0.52337644549944951],
    [-4.17871859155190428, -0.32768282076106237, 0.47662355450055044],
    [0.50287263494578682, -2.57192694985560522, 0.59603920482822492]
])
# This linear combinations are used in the algorithm.
TI_REAL = TI[0]
TI_COMPLEX = TI[1] + 1j * TI[2]

NEWTON_MAXITER = 7  # Maximum number of Newton iterations.
MAX_FACTOR = 8  # Maximum allowed increase in a step size.
MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.


def create_spline(x, y, Z):
    """Create a cubic spline given values at 4 points on each interval.

    Parameters
    ----------
    x : ndarray, shape (n_points,)
        Values of the independent variable.
    y : ndarray, shape (n_points, n)
        Values of the dependent variable at `x`.
    Z : ndarray, shape (n_points - 1, 3, n)
        Values to add to `y` to compute the solution at 3 other points on
        each interval (see `C` constant vector at the beginning of this file).

    Returns
    -------
    sol : PPoly
        Constructed spline as a PPoly instance.
    """
    from scipy.interpolate import PPoly

    if x[-1] < x[0]:
        reverse = True
        x = x[::-1]
        y = y[::-1]
        Z = Z[::-1]

        z0 = Z[:, 1] - Z[:, 2]
        z1 = Z[:, 0] - Z[:, 2]
        z2 = -Z[:, 2]
    else:
        reverse = False
        z0 = Z[:, 0]
        z1 = Z[:, 1]
        z2 = Z[:, 2]

    h = np.diff(x)[:, None]
    n_points, n = y.shape
    c = np.empty((4, n_points - 1, n))

    if reverse:
        c[0] = ((-10 + 15*S6) * z0 - (10 + 15*S6) * z1 + 30 * z2) / (3 * h**3)
        c[1] = ((7 - 23*S6) * z0 + (7 + 23*S6) * z1 - 36 * z2) / (3 * h**2)
        c[2] = ((1 + 8*S6/3) * z0 + (1 - 8*S6/3) * z1 + 3 * z2) / h
    else:
        c[0] = ((10 + 15*S6) * z0 + (10 - 15*S6) * z1 + 10 * z2) / (3 * h**3)
        c[1] = -((23 + 22*S6) * z0 + (23 - 22*S6) * z1 + 8 * z2) / (3 * h**2)
        c[2] = ((13 + 7*S6) * z0 + (13 - 7*S6) * z1 + z2) / (3 * h)
    c[3] = y[:-1]

    c = np.rollaxis(c, 2)
    return PPoly(c, x, extrapolate=True, axis=1)


def create_spline_one_step(x, x_new, y, Z):
    """Create a cubic spline for a single step.

    Parameters
    ----------
    x, x_new : float
        Previous and new values of the independed variable.
    y : ndarray, shape (n,)
        Previous value of the dependent variable.
    Z : ndarray, shape (3, n)
        Values to add to `y` to compute the solution at 3 other points on
        the interval (see `C` constant vector at the beginning of this file).

    Returns
    -------
    sol : PPoly
        Constructed spline as a PPoly instance.
    """

    from scipy.interpolate import PPoly

    if x_new < x:
        reverse = True
        x, x_new = x_new, x
        y = y + Z[2]
        z0 = Z[1] - Z[2]
        z1 = Z[0] - Z[2]
        z2 = -Z[2]
    else:
        reverse = False
        z0 = Z[0]
        z1 = Z[1]
        z2 = Z[2]

    h = x_new - x
    n = y.shape[0]
    c = np.empty((4, 1, n))

    if reverse:
        c[0] = ((-10 + 15*S6) * z0 - (10 + 15*S6) * z1 + 30 * z2) / (3 * h**3)
        c[1] = ((7 - 23*S6) * z0 + (7 + 23*S6) * z1 - 36 * z2) / (3 * h**2)
        c[2] = ((1 + 8*S6/3) * z0 + (1 - 8*S6/3) * z1 + 3 * z2) / h
    else:
        c[0] = ((10 + 15*S6) * z0 + (10 - 15*S6) * z1 + 10 * z2) / (3 * h**3)
        c[1] = -((23 + 22*S6) * z0 + (23 - 22*S6) * z1 + 8 * z2) / (3 * h**2)
        c[2] = ((13 + 7*S6) * z0 + (13 - 7*S6) * z1 + z2) / (3 * h)
    c[3] = y

    return PPoly(c, [x, x_new], extrapolate=True)


def solve_collocation_system(fun, x, y, h, J, Z0, scale, tol, LU_real,
                             LU_complex):
    """Solve the collocation system.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    x : float
        Current value of the independent variable.
    y : ndarray, shape (n,)
        Current value of the dependent variable.
    h : float
        Step to try.
    J : ndarray, shape (n, n)
        Jacobian of `fun` with respect to `y`.
    Z0 : ndarray, shape (3, n)
        Initial guess for the solution.
    scale : float
        Problem tolerance scale, i.e. ``rtol * abs(y) + atol``.
    tol : float
        Tolerance to which solve the system.
    LU_real, LU_complex : tuple or None
        LU decompositions of the system Jacobian. If None, they will be
        compute inside the function.

    Returns
    -------
    status : int
        Status of the solution:

            * 0: Iterations converged.
            * 1: Iterations didn't converge.

        Potentially to be extended with the status when the system Jacobian
        is singular.
    n_iter : int
        Number of completed iterations.
    Z : ndarray, shape (3, n)
        Found solution.
    f_new : ndarray, shape (3, n)
        Value of `fun(x + h, y(x + h))`.
    theta : float
        The rate of convergence.
    LU_real, LU_complex : tuple
        Computed LU decompositions.
    """
    n = y.shape[0]
    I = np.identity(n)
    M_real = MU_REAL / h
    M_complex = MU_COMPLEX / h

    if LU_real is None or LU_complex is None:
        LU_real = lu_factor(M_real * I - J)
        LU_complex = lu_factor(M_complex * I - J)

    W = TI.dot(Z0)
    Z = Z0

    F = np.empty((3, n))
    ch = h * C

    dW_norm_old = None
    dW = np.empty_like(W)
    status = 0
    for k in range(NEWTON_MAXITER):
        for i in range(3):
            F[i] = fun(x + ch[i], y + Z[i])

        f_real = F.T.dot(TI_REAL) - M_real * W[0]
        f_complex = F.T.dot(TI_COMPLEX) - M_complex * (W[1] + 1j * W[2])

        dW_real = lu_solve(LU_real, f_real, overwrite_b=True)
        dW_complex = lu_solve(LU_complex, f_complex, overwrite_b=True)

        dW[0] = dW_real
        dW[1] = dW_complex.real
        dW[2] = dW_complex.imag

        dW_norm = norm(dW / scale)
        if dW_norm_old is not None:
            theta = dW_norm / dW_norm_old
            eta = theta / (1 - theta)
        else:
            theta = 0
            eta = 1

        if (theta > 1 or
                theta ** (NEWTON_MAXITER - k) / (1 - theta) * dW_norm > tol):
            status = 1
            break

        W += dW
        Z = T.dot(W)

        if eta * dW_norm < tol:
            break

        dW_norm_old = dW_norm

    return status, k + 1, Z, F[-1], theta, LU_real, LU_complex


def predict_factor(h_abs, h_abs_old, error_norm, error_norm_old):
    """Predict by which factor to increase the step size.

    The algorithm is described in [1]_.

    Parameters
    ----------
    h_abs, h_abs_old : float
        Current and previous values of the step size, `h_abs_old` can be None
        (see Notes).
    error_norm, error_norm_old : float
        Current and previous values of the error norm, `error_norm_old` can
        be None (see Notes).

    Returns
    -------
    factor : float
        Predicted factor.

    Notes
    -----
    If `h_abs_old` and `error_norm_old` are both not None then a two-step
    algorithm is used, otherwise a one-step algorithm is used.

    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations II: Stiff and Differential-Algebraic Problems", Sec. IV.8.
    """
    with np.errstate(divide='ignore'):
        if error_norm_old is None or h_abs_old is None:
            multiplier = 1
        else:
            multiplier = h_abs / h_abs_old * (error_norm_old /
                                              error_norm) ** 0.25

        factor = min(1, multiplier) * error_norm ** -0.25

    return factor


def radau(fun, jac, a, b, ya, fa, Ja, rtol, atol, max_step,
          events, is_terminal, direction):
    """
    Integrate an ODE by an implicit Runge-Kutta method of Radau IIA family.
    """
    s = np.sign(b - a)
    h_abs = select_initial_step(fun, a, b, ya, fa, 5, rtol, atol)

    newton_tol = max(10 * EPS / rtol, min(0.03, rtol ** 0.5))

    x = a
    y = ya
    f = fa
    J = Ja
    current_jac = True

    xs = [x]
    ys = [y]
    fs = [f]
    Zs = []

    if events is not None:
        g = [event(x, y) for event in events]
        x_events = [[] for _ in range(len(events))]
    else:
        x_events = None

    LU_real = None
    LU_complex = None

    sol = None
    rejected = True

    h_abs_old = None
    error_norm_old = None

    status = None
    while status is None:
        if h_abs > max_step:
            h_abs = max_step
            h_abs_old = None
            error_norm_old = None

        d = abs(b - x)
        if h_abs > d:
            status = 0
            h_abs = d
            x_new = b
            h = h_abs * s
            h_abs_old = None
            error_norm_old = None
        else:
            h = h_abs * s
            x_new = x + h
            if x_new == x:  # h less than spacing between numbers.
                status = -1

        if sol is None:
            Z0 = np.zeros((3, y.shape[0]))
        else:
            Z0 = sol(x + h * C) - y

        scale = atol + np.abs(y) * rtol
        newton_status, n_iter, Z, f_new, theta, LU_real, LU_complex = \
            solve_collocation_system(fun, x, y, h, J, Z0, scale,
                                     newton_tol, LU_real, LU_complex)
        if newton_status == 1:
            status = None
            rejected = True
            if not current_jac:
                J = jac(x, y)
                h_abs *= 0.75
            else:
                h_abs *= 0.5
            LU_real = None
            LU_complex = None
            continue

        y_new = y + Z[-1]

        ZE = Z.T.dot(E) / h
        error = lu_solve(LU_real, f + ZE)
        scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
        error_norm = norm(error / scale)
        safety = 0.9 * (2 * NEWTON_MAXITER + 1) / (2 * NEWTON_MAXITER + n_iter)

        if rejected and error_norm > 1:
            error = lu_solve(LU_real, fun(x, y + error) + ZE)
            error_norm = norm(error / scale)

        if error_norm > 1:
            factor = predict_factor(h_abs, h_abs_old,
                                    error_norm, error_norm_old)
            h_abs *= max(MIN_FACTOR, safety * factor)

            status = None
            rejected = True
            LU_real = None
            LU_complex = None
            continue

        sol = create_spline_one_step(x, x_new, y, Z)
        if events is not None:
            g_new = [event(x_new, y_new) for event in events]
            active_events = get_active_events(g, g_new, direction)
            g = g_new
            if active_events.size > 0:
                root_indices, roots, terminate = handle_events(
                    sol, events, active_events, is_terminal, x, x_new)

                for e, xe in zip(root_indices, roots):
                    x_events[e].append(xe)

                if terminate:
                    status = 1
                    x_new = roots[-1]
                    y_new = sol(x_new)
                    f_new = fun(x_new, y_new)
                    Z = sol(x + (x_new - x) * C) - y

        recompute_jac = jac is not None and n_iter > 1 and theta > 1e-3

        factor = predict_factor(h_abs, h_abs_old, error_norm, error_norm_old)
        factor = min(MAX_FACTOR, safety * factor)

        if not recompute_jac and factor < 1.2:
            factor = 1
        else:
            LU_real = None
            LU_complex = None

        if recompute_jac:
            J = jac(x, y)
            current_jac = True
        elif jac is not None:
            current_jac = False

        h_abs_old = h_abs
        h_abs *= factor
        error_norm_old = error_norm

        x = x_new
        y = y_new
        f = f_new

        xs.append(x)
        ys.append(y)
        fs.append(f)
        Zs.append(Z)

        rejected = False

    xs = np.asarray(xs)
    ys = np.asarray(ys)
    fs = np.asarray(fs)
    Zs = np.asarray(Zs)

    sol = create_spline(xs, ys, Zs)

    if x_events:
        x_events = [np.asarray(xe) for xe in x_events]
        if len(x_events) == 1:
            x_events = x_events[0]

    return status, sol, xs, ys.T, fs.T, x_events
