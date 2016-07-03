"""Explicit Runge-Kutta methods."""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.interpolate import PPoly
from .common import select_initial_step, get_active_events, handle_events, norm

# Algorithm parameters.

# Multiply steps computed from asymptotic behaviour of errors by this.
SAFETY = 0.9

MAX_FACTOR = 5  # Maximum allowed increase in a step size.
MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.

# Butcher tables. See `rk_step` for explanation.

# Bogacki–Shampine scheme.
C23 = np.array([1/2, 3/4])
A23 = [np.array([1/2]),
       np.array([0, 3/4])]
B23 = np.array([2/9, 1/3, 4/9])
# Coefficients for estimation errors. The difference between B's for lower
# and higher order accuracy methods.
E23 = np.array([5/72, -1/12, -1/9, 1/8])

# Dormand–Prince scheme.
C45 = np.array([1/5, 3/10, 4/5, 8/9, 1])
A45 = [np.array([1/5]),
       np.array([3/40, 9/40]),
       np.array([44/45, -56/15, 32/9]),
       np.array([19372/6561, -25360/2187, 64448/6561, -212/729]),
       np.array([9017/3168, -355/33, 46732/5247, 49/176, -5103/18656])]
B45 = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84])
E45 = np.array([-71/57600, 0, 71/16695, -71/1920, 17253/339200, -22/525, 1/40])

# Coefficients to compute y(x + 0.5 * h) from RK stages with a 4-rd order
# accuracy. Then it can be used for quartic interpolation with a 4-rd order
# accuracy.
M45 = np.array([613/3072, 0, 125/159, -125/1536, 8019/54272, -11/96, 1/16])


def prepare_method(method, n, s):
    """Choose appropriate matrices for a RK method.

    See `rk_step` for the explanation of returned matrices.
    """
    if method == 'RK45':
        A = A45
        B = B45
        C = C45
        E = E45
        M = M45
        K = np.empty((7, n))
        order = 5
    elif method == 'RK23':
        A = A23
        B = B23
        C = C23
        E = E23
        M = None
        order = 3
        K = np.empty((4, n))
    else:
        raise ValueError("`method` must be 'RK45' or 'RK23'.")

    return A, B, C, E, M, K, order


def rk_step(fun, x, y, f, h, A, B, C, E, K):
    """Perform a single Runge-Kutta step.

    This function computes a prediction of an explicit Runge-Kutta method and
    also estimates the error of a less accurate method.

    Notation for Butcher tableau is as in [1]_.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    x : float
        Current value of the independent variable.
    y : ndarray, shape (n,)
        Current value of the solution.
    f : ndarray, shape (n,)
        Current value of the derivative of the solution, i.e. ``fun(x, y)``.
    h : float, shape (n,)
        Step for x to use.
    A : list of ndarray, length n_stages - 1
        Coefficients for combining previous RK stages for computing the next
        stage. For explicit methods the coefficients above the main diagonal
        are zeros, so they are stored as a list of arrays of increasing
        lengths. The first stage is always just `f`, thus no coefficients are
        required.
    B : ndarray, shape (n_stages,)
        Coefficients for combining RK stages for computing the final
        prediction.
    C : ndarray, shape (n_stages - 1,)
        Coefficients for incrementing x for computing RK stages. The value for
        the first stage is always zero, thus it is not stored.
    E : ndarray, shape (n_stages + 1,)
        Coefficients for estimating the error of a less accurate method. They
        are computed as the difference between b's in an extended tableau.
    K : ndarray, shape (n_stages + 1, n)
        Storage array for putting RK stages here. Stages are stored in rows.

    Returns
    -------
    y_new : ndarray, shape (n,)
        Solution at x + h computed with a higher accuracy.
    f_new : ndarray, shape (n,)
        Derivative ``fun(x + h, y_new)``.
    error : ndarray, shape (n,)
        Error estimate.

    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations I: Nonstiff Problems", Sec. II.4.
    """
    K[0] = f
    for s, (a, c) in enumerate(zip(A, C)):
        dy = np.dot(K[:s + 1].T, a) * h
        K[s + 1] = fun(x + c * h, y + dy)

    y_new = y + h * np.dot(K[:-1].T, B)
    f_new = fun(x + h, y_new)

    K[-1] = f_new
    error = np.dot(K.T, E) * h

    return y_new, f_new, error


def create_spline(x, y, f, ym):
    """Create a cubic or quartic spline given values and derivatives.

    Parameters
    ----------
    x : ndarray, shape (n_points,)
        Values of the independent variable.
    y : ndarray, shape (n_points, n)
        Values of the dependent variable at `x`.
    f : ndarray, shape (n_points, n)
        Values of the derivatives of `y` evaluated at `x`.
    ym : ndarray with shape (n_points, n) or None
        Values of the dependent variables at middle points between values
        of `x`. If None, a cubic spline will be constructed, and a quartic
        spline otherwise.

    Returns
    -------
    sol : PPoly
        Constructed spline as a PPoly instance.
    """
    from scipy.interpolate import PPoly

    if x[-1] < x[0]:
        x = x[::-1]
        y = y[::-1]
        if ym is not None:
            ym = ym[::-1]
        f = f[::-1]

    h = np.diff(x)

    y0 = y[:-1]
    y1 = y[1:]
    f0 = f[:-1]
    f1 = f[1:]

    n_points, n = y.shape
    h = h[:, None]
    if ym is None:
        c = np.empty((4, n_points - 1, n))
        slope = (y1 - y0) / h
        t = (f0 + f1 - 2 * slope) / h
        c[0] = t / h
        c[1] = (slope - f0) / h - t
        c[2] = f0
        c[3] = y0
    else:
        c = np.empty((5, n_points - 1, n))
        c[0] = (-8 * y0 - 8 * y1 + 16 * ym) / h**4 + (- 2 * f0 + 2 * f1) / h**3
        c[1] = (18 * y0 + 14 * y1 - 32 * ym) / h**3 + (5 * f0 - 3 * f1) / h**2
        c[2] = (-11 * y0 - 5 * y1 + 16 * ym) / h**2 + (-4 * f0 + f1) / h
        c[3] = f0
        c[4] = y0

    c = np.rollaxis(c, 2)
    return PPoly(c, x, extrapolate=True, axis=1)


def create_spline_one_step(x, x_new, y, y_new, f, f_new, ym):
    """Create a spline for a single step.

    Parameters
    ----------
    x, x_new : float
        Previous and new values of the independed variable.
    y, y_new : float
        Previous and new values of the dependent variable.
    f, f_new : float
        Previous and new values of the derivative of the dependent variable.
    ym : float or None
        Value of the dependent variable at the middle point between `x` and
        `x_new`. If provided the quartic spline is constructed, if None
        the cubic spline is constructed.

    Returns
    -------
    sol : PPoly
        Constructed spline as a PPoly instance.
    """
    if x_new < x:
        x0, x1 = x_new, x
        y0, y1 = y_new, y
        f0, f1 = f_new, f
    else:
        x0, x1 = x, x_new
        y0, y1 = y, y_new
        f0, f1 = f, f_new

    h = x1 - x0
    n = y.shape[0]
    if ym is None:
        c = np.empty((4, 1, n), dtype=y.dtype)
        slope = (y1 - y0) / h
        t = (f0 + f1 - 2 * slope) / h
        c[0] = t / h
        c[1] = (slope - f0) / h - t
        c[2] = f0
        c[3] = y0
    else:
        c = np.empty((5, 1, n), dtype=y.dtype)
        c[0] = (-8 * y0 - 8 * y1 + 16 * ym) / h**4 + (- 2 * f0 + 2 * f1) / h**3
        c[1] = (18 * y0 + 14 * y1 - 32 * ym) / h**3 + (5 * f0 - 3 * f1) / h**2
        c[2] = (-11 * y0 - 5 * y1 + 16 * ym) / h**2 + (-4 * f0 + f1) / h
        c[3] = f0
        c[4] = y0

    c = np.rollaxis(c, 2)
    return PPoly(c, [x0, x1], extrapolate=True, axis=1)


def rk(fun, a, b, ya, fa, rtol, atol, method, events, direction, is_terminal):
    """Integrate an ODE by Runge-Kutta method."""
    max_step = 0.1 * np.abs(b - a)
    s = np.sign(b - a)

    A, B, C, E, M, K, order = prepare_method(method, ya.shape[0], s)

    h_abs = select_initial_step(fun, a, b, ya, fa, order, rtol, atol)

    x = a
    y = ya
    f = fa

    ys = [y]
    xs = [x]
    fs = [f]
    if order == 3:
        yms = None
    else:
        yms = []

    if events is not None:
        g = [event(x, y) for event in events]
        x_events = [[] for _ in range(len(events))]
    else:
        x_events = None

    status = None
    while status is None:
        h_abs = min(h_abs, max_step)

        d = abs(b - x)
        if h_abs > d:
            status = 1
            h_abs = d
            x_new = b
            h = h_abs * s
        else:
            h = h_abs * s
            x_new = x + h
            if x_new == x:  # h less than spacing between numbers.
                status = 0

        y_new, f_new, error = rk_step(fun, x, y, f, h, A, B, C, E, K)
        scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
        error_norm = norm(error / scale)

        if error_norm > 1:
            h_abs *= max(MIN_FACTOR, SAFETY * error_norm**(-1/order))
            status = None
            continue

        if M is not None:
            ym = y + 0.5 * h * np.dot(K.T, M)
        else:
            ym = None

        if events is not None:
            g_new = [event(x_new, y_new) for event in events]
            active_events = get_active_events(g, g_new, direction)
            g = g_new
            if active_events.size > 0:
                sol = create_spline_one_step(x, x_new, y, y_new, f, f_new, ym)
                root_indices, roots, terminate = handle_events(
                    sol, events, active_events, is_terminal, x, x_new)

                for e, xe in zip(root_indices, roots):
                    x_events[e].append(xe)

                if terminate:
                    status = 2
                    x_new = roots[-1]
                    y_new = sol(x_new)
                    if ym is not None:
                        ym = sol(0.5 * (x + x_new))
                    f_new = fun(x_new, y_new)

        x = x_new
        y = y_new
        f = f_new
        ys.append(y)
        xs.append(x)
        fs.append(f)
        if ym is not None:
            yms.append(ym)

        with np.errstate(divide='ignore'):
            h_abs *= min(MAX_FACTOR, SAFETY * error_norm**(-1/order))

    xs = np.asarray(xs)
    ys = np.asarray(ys)
    fs = np.asarray(fs)
    if yms is not None:
        yms = np.asarray(yms)

    sol = create_spline(xs, ys, fs, yms)

    if x_events:
        x_events = [np.asarray(xe) for xe in x_events]
        if len(x_events) == 1:
            x_events = x_events[0]

    return status, sol, xs, ys.T, fs.T, x_events
