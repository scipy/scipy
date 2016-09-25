from __future__ import division, print_function, absolute_import
import numpy as np
from .bdf import BDF
from .radau import Radau
from .rk import RK23, RK45
from scipy.optimize import OptimizeResult
from .common import EPS, OdeSolution


METHODS = ['RK23', 'RK45', 'Radau', 'BDF']


class OdeResult(OptimizeResult):
    pass


def prepare_events(events):
    if callable(events):
        events = (events,)

    if events is not None:
        is_terminal = np.empty(len(events), dtype=bool)
        direction = np.empty(len(events))
        for i, event in enumerate(events):
            try:
                is_terminal[i] = event.terminal
            except AttributeError:
                is_terminal[i] = False

            try:
                direction[i] = event.direction
            except AttributeError:
                direction[i] = 0
    else:
        is_terminal = None
        direction = None

    return events, is_terminal, direction


def solve_event_equation(event, sol, t_old, t):
    """Solve an equation corresponding to an ODE event.

    The equation is ``event(t, y(t)) = 0``, here ``y(t)`` is known from an
    ODE solver using some sort of interpolation. It is solved by
    `scipy.optimize.brentq` with xtol=atol=4*EPS.

    Parameters
    ----------
    event : callable
        Function ``event(x, y)``.
    sol : callable
        Computed solution ``y(x)``. It should be defined only between `x` and
        `x_new`.
    t_old, t : float
        Previous and new values of the independed variable, it will be used as
        a bracketing interval.

    Returns
    -------
    root : float
        Found solution.
    """
    from scipy.optimize import brentq
    return brentq(lambda t: event(t, sol(t)), t_old, t, xtol=4 * EPS)


def handle_events(sol, events, active_events, is_terminal, t_old, t):
    """Helper function to handle events.

    Parameters
    ----------
    sol : DenseOutput
        Function ``sol(x)`` which evaluates an ODE solution.
    events : list of callables, length n_events
        Event functions.
    active_events : ndarray
        Indices of events which occurred
    is_terminal : ndarray, shape (n_events,)
        Which events are terminate.
    t_old, t : float
        Previous and new values of time.

    Returns
    -------
    root_indices : ndarray
        Indices of events which take zero before a possible termination.
    roots : ndarray
        Values of x at which events take zero values.
    terminate : bool
        Whether a termination event occurred.
    t_old, t : float
        Previous and new values of the independed variable, it will be used as
        a bracketing interval.
    """
    roots = []
    for event_index in active_events:
        roots.append(solve_event_equation(events[event_index], sol, t_old, t))

    roots = np.asarray(roots)

    if np.any(is_terminal[active_events]):
        if t > t_old:
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


def find_active_events(g, g_new, direction):
    """Find which event occurred during an integration step.

    Parameters
    ----------
    g, g_new : array_like, shape (n_events,)
        Values of event functions at a current and next points.
    direction : ndarray, shape (n_events,)
        Event "direction" according to the definition in `solve_ivp`.

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


def solve_ivp(fun, t_span, y0, rtol=1e-3, atol=1e-6, method='RK45',
              dense_output=False, max_step=np.inf, jac=None, events=None):
    """Solve an initial value problem for a system of ODEs.

    This function numerically integrates a system of ODEs given an initial
    value::

        dy / dt = f(t, y)
        y(t0) = y0

    Here t is a 1-dimensional independent variable (time), y(t) is a
    n-dimensional vector-valued function (state) and y0 is an initial state.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is ``fun(t, y)``.
        Here ``t`` is a scalar, and ``y`` is ndarray with shape (n,). It
        must return an array_like with shape (n,).
    t_span : 2-tuple of floats
        Interval of integration (t0, tf). The solver starts with t=t0 and
        integrates until it reaches t=tf.
    y0 : array_like, shape (n,)
        Initial state.
    rtol, atol : float and array_like, optional
        Relative and absolute tolerances. The solver keeps the local error
        estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
        relative accuracy (number of correct digits). But if a component of `y`
        is approximately below `atol` then the error only needs to fall within
        the same `atol` threshold, and the number of correct digits is not
        guaranteed. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    method : string, optional
        Integration method to use:

            * 'RK45' (default): Explicit Runge-Kutta method of order 5(4) [1]_.
              The error is controlled assuming 4th order accuracy, but steps
              are taken using a 5th oder accurate formula (local extrapolation
              is done). A quartic interpolation polynomial is used for the
              dense output [2]_.
            * 'RK23': Explicit Runge-Kutta method of order 3(2) [3]_. The error
              is controlled assuming 2nd order accuracy, but steps are taken
              using a 3rd oder accurate formula (local extrapolation is done).
              A cubic Hermit polynomial is used for the dense output.
            * 'Radau': Implicit Runge-Kutta method of Radau IIA family of
              order 5 [4]_. The error is controlled for a 3rd order accurate
              embedded formula. A cubic polynomial which satisfies the
              collocation conditions is used for the dense output.
            * 'BDF': Implicit multi-step variable order (1 to 5) method based
              on a Backward Differentiation Formulas for the derivative
              approximation [5]_. An implementation approach follows the one
              described in [6]_. A quasi-constant step scheme is used
              and accuracy enhancement using NDF modification is also
              implemented.

        You should use 'RK45' or 'RK23' methods for non-stiff problems and
        'Radau' or 'BDF' for stiff problems [7]_. If not sure, first try to run
        'RK45' and if it does unusual many iterations or diverges then your
        problem is likely to be stiff and you should use 'Radau' or 'BDF'.
    dense_output : bool, optional
        Whether to compute a continuous solution. Default is False.
    max_step : float, optional
        Maximum allowed step size. Default is np.inf, i.e. step is not
        bounded and determined solely by the solver.
    jac : array_like, callable or None, optional
        Jacobian matrix of the right-hand side of the system with respect to
        `y`, required only by 'Radau' and 'BDF' methods. The Jacobian matrix
        has shape (n, n) and its element (i, j) is equal to ``d f_i / d y_j``.
        There are 3 ways to define the Jacobian:

            * If array_like, then the Jacobian is assumed to be constant.
            * If callable, then the Jacobian is assumed to depend on both
              t and y, and will be called as ``jac(t, y)`` as necessary.
            * If None (default), then the Jacobian will be approximated by
              finite differences.

        It is generally recommended to provided the Jacobian rather than
        relying on finite difference approximation.
    events : callable, list of callables or None, optional
        Events to track. Events are defined by functions which take
        a zero value at a point of an event. Each function must have a
        signature ``event(t, y)`` and return float, the solver will find an
        accurate value of ``t`` at which ``event(t, y(t)) = 0`` using a root
        finding algorithm. Additionally each ``event`` function might have
        attributes:

            * terminal: bool, whether to terminate integration if this
              event occurs. Implicitly False if not assigned.
            * direction: float, direction of crossing a zero. If `direction`
              is positive then `event` must go from negative to positive, and
              vice-versa if `direction` is negative. If 0, then either way will
              count. Implicitly 0 if not assigned.

        You can assign attributes like ``event.terminal = True`` to any
        function in Python. If None (default), events won't be tracked.

    Returns
    -------
    Bunch object with the following fields defined:
    t : ndarray, shape (n_points,)
        Time points at which the solver made steps.
    y : ndarray, shape (n, n_points)
        Solution values at `t`.
    sol : `OdeSolution` or None
        Found solution as `OdeSolution` instance, None if `dense_output` was
        set to False.
    t_events : list of ndarray or None
        Contains arrays with times at each a corresponding event was detected,
        the length of the list equals to the number of events. None if `events`
        was None.
    status : int
        Reason for algorithm termination:

            * -1: Integration step failed.
            * 0: The solver successfully reached the interval end.
            * 1: A termination event occurred.

    message : string
        Verbal description of the termination reason.
    success : bool
        True if the solver reached the interval end or a termination event
        occurred (``status >= 0``).

    References
    ----------
    .. [1] J. R. Dormand, P. J. Prince, "A family of embedded Runge-Kutta
           formulae", Journal of Computational and Applied Mathematics, Vol. 6,
           No. 1, pp. 19-26, 1980.
    .. [2] L. W. Shampine, "Some Practical Runge-Kutta Formulas", Mathematics
           of Computation,, Vol. 46, No. 173, pp. 135-150, 1986.
    .. [3] P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
           Appl. Math. Lett. Vol. 2, No. 4. pp. 321-325, 1989.
    .. [4] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations I:
           Nonstiff Problems", Sec. III.1.
    .. [5] `Backward Differentiation Formula
            <https://en.wikipedia.org/wiki/Backward_differentiation_formula>`_
            on Wikipedia.
    .. [6] L. F. Shampine, M. W. Reichelt, "THE MATLAB ODE SUITE", SIAM J. SCI.
           COMPUTE., Vol. 18, No. 1, pp. 1-22, January 1997.
    .. [7] `Stiff equation <https://en.wikipedia.org/wiki/Stiff_equation>`_ on
           Wikipedia.
    """
    if method not in METHODS:
        raise ValueError("`method` must be one of {}.".format(METHODS))

    t0, tf = float(t_span[0]), float(t_span[1])
    if t0 == tf:
        raise ValueError("Initial and final times must be distinct.")

    y0 = np.atleast_1d(y0)
    if y0.ndim != 1:
        raise ValueError("`y0` must be 1-dimensional.")

    if method == 'RK23':
        solver = RK23(fun, t0, y0, tf, rtol, atol)
    elif method == 'RK45':
        solver = RK45(fun, t0, y0, tf, rtol, atol)
    elif method == 'Radau':
        solver = Radau(fun, t0, y0, tf, rtol, atol, jac)
    elif method == 'BDF':
        solver = BDF(fun, t0, y0, tf, rtol, atol, jac)

    ts = [t0]
    ys = [y0]
    interpolants = []

    events, is_terminal, event_dir = prepare_events(events)

    if events is not None:
        g = [event(t0, y0) for event in events]
        t_events = [[] for _ in range(len(events))]
    else:
        t_events = None

    while solver.status != 'finished':
        status, message = solver.step(max_step)
        if status == 'failed':
            break

        t = solver.t
        y = solver.y

        ts.append(t)
        ys.append(y)

        if dense_output:
            sol = solver.dense_output()
            interpolants.append(sol)
        else:
            sol = None

        if events is not None:
            g_new = [event(t, y) for event in events]
            active_events = find_active_events(g, g_new, event_dir)
            if active_events.size > 0:
                if sol is None:
                    sol = solver.dense_output()
                root_indices, roots, terminate = handle_events(
                    sol, events, active_events, is_terminal, ts[-2], ts[-1])

                for e, te in zip(root_indices, roots):
                    t_events[e].append(te)

                if terminate:
                    status = 1
                    t_term = roots[-1]
                    ts[-1] = t_term
                    ys[-1] = sol(t_term)
                    break
            g = g_new

    if status == 'finished':
        status = 0
        message = "The solver successfully reached the interval end."
    elif status == 'running':
        status = 1
        message = "A termination event occurred."
    elif status == 'failed':
        status = -1

    if t_events is not None:
        t_events = [np.asarray(xe) for xe in t_events]

    ts = np.array(ts)
    ys = np.array(ys).T

    if dense_output:
        sol = OdeSolution(ts, interpolants)
    else:
        sol = None

    return OdeResult(t=ts, y=ys, sol=sol, t_events=t_events,
                     status=status, message=message, success=status >= 0)
