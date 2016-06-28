"""Common functions for ODE solvers."""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.optimize import brentq, OptimizeResult


EPS = np.finfo(float).eps


class ODEResult(OptimizeResult):
    pass


def norm(x):
    """Compute RMS norm."""
    return np.linalg.norm(x) / x.size ** 0.5


def select_initial_step(fun, a, b, ya, fa, order, rtol, atol):
    """Empirically select a good initial step.

    The algorithm is described in [1]_.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    a : float
        Initial value of the independent variable.
    b : float
        Final value value of the independent variable.
    ya : ndarray, shape (n,)
        Initial value of the dependent variable.
    fa : ndarray, shape (n,)
        Initial value of the derivative, i. e. ``fun(x0, y0)``.
    order : float
        Method order.
    rtol : float
        Desired relative tolerance.
    atol : float
        Desired absolute tolerance.

    Returns
    -------
    h_abs : float
        Absolute value of the suggested initial step.

    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations I: Nonstiff Problems", Sec. II.4.
    """
    scale = atol + np.abs(ya) * rtol
    d0 = norm(ya / scale)
    d1 = norm(fa / scale)
    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1

    s = np.sign(b - a)
    y1 = ya + h0 * s * fa
    f1 = fun(a + h0 * s, y1)
    d2 = norm((f1 - fa) / scale) / h0

    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1 / order)

    return min(100 * h0, h1)


def get_active_events(g, g_new, direction):
    """Find which event occurred during an integration step.

    Parameters
    ----------
    g, g_new : array_like, shape (n_events,)
        Values of event functions at a current and next points.
    direction : ndarray, shape (n_events,)
        Event "direction" according to definition in `solve_ivp`.

    Returns
    -------
    active_events : ndarray
        Indices of events which occurred during the step.
    """
    g, g_new = np.asarray(g), np.asarray(g_new)
    up = (g <= 0) & (g_new >= 0)
    down = (g >= 0) & (g_new <= 0)
    either = up | down
    mask = (up & (direction > 0) |
            down & (direction < 0) |
            either & (direction == 0))

    return np.nonzero(mask)[0]


def handle_events(sol, events, active_events, is_terminal, x, x_new):
    """Helper function to handle events.

    Parameters
    ----------
    sol : callable
        Function ``sol(x)`` which evaluates an ODE solution.
    events : list of callables, length n_events
        Event functions.
    active_events : ndarray
        Indices of events which occurred
    is_terminal : ndarray, shape (n_events,)
        Which events are terminate.
    x, x_new : float
        Previous and new values of the independed variable, it will be used as
        a bracketing interval.

    Returns
    -------
    root_indices : ndarray
        Indices of events which take zero before a possible termination.
    roots : ndarray
        Values of x at which events take zero values.
    terminate : bool
        Whether a termination event occurred.
    """
    roots = []
    for event_index in active_events:
        roots.append(solve_event_equation(events[event_index], sol, x, x_new))

    roots = np.asarray(roots)

    if np.any(is_terminal[active_events]):
        if x_new > x:
            order = np.argsort(roots)
        else:
            order = np.argsort(-roots)
        active_events = active_events[order]
        roots = roots[order]
        t = np.nonzero(is_terminal[active_events])[0][0]
        active_events = active_events[:t + 1]
        roots = roots[:t + 1]
        terminate = True
    else:
        terminate = False

    return active_events, roots, terminate


def solve_event_equation(event, sol, x, x_new):
    """Solve an equation corresponding to an ODE event.

    The equation is ``event(x, y(x)) = 0``, here ``y(x)`` is known from an
    ODE solver using some sort of interpolation. It is solved by
    `scipy.optimize.brentq` with xtol=atol=4*EPS.

    Parameters
    ----------
    event : callable
        Function ``event(x, y)``.
    sol : callable
        Computed solution ``y(x)``. It should be defined only between `x` and
        `x_new`.
    x, x_new : float
        Previous and new values of the independed variable, it will be used as
        a bracketing interval.

    Returns
    -------
    root : float
        Found solution.
    """
    return brentq(lambda t: event(t, sol(t)), x, x_new, xtol=4 * EPS)
