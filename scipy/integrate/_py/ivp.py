"""Generic interface for initial value problem solvers."""
from __future__ import division, print_function, absolute_import

from warnings import warn
import numpy as np
from .common import select_initial_step, EPS, ODEResult
from .rk import rk


METHODS = ['RK23', 'RK45']


TERMINATION_MESSAGES = {
    0: "The solver failed to reach the interval end or a termination event.",
    1: "The solver successfully reached the interval end.",
    2: "A termination event occurred."
}


def validate_tol(rtol, atol, n):
    """Validate tolerance values."""
    if rtol < 100 * EPS:
        warn("`rtol` is too low, setting to {}".format(100 * EPS))
        rtol = 100 * EPS

    atol = np.asarray(atol)
    if atol.ndim > 0 and atol.shape != (n,):
        raise ValueError("`atol` has wrong shape.")

    if np.any(atol < 0):
        raise ValueError("`atol` must be positive.")

    return rtol, atol


def solve_ivp(fun, x_span, ya, rtol=1e-3, atol=1e-6, method='RK45',
              events=None):
    """Solve an initial value problem for a system of ODEs.

    This function numerically integrates a system of ODEs given an initial
    value::

        dy / dx = f(x, y)
        y(a) = ya

    Here x is a 1-dimensional independent variable, y(x) is a n-dimensional
    vector-valued function and ya is a n-dimensional vector with initial
    values.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is ``fun(x, y)``.
        Here ``x`` is a scalar, and ``y`` is ndarray with shape (n,). It
        must return an array_like with shape (n,).
    x_span : 2-tuple of floats
        Interval of integration (a, b). The solver starts with x=a and
        integrates until it reaches x=b.
    ya : array_like, shape (n,)
        Initial values for y.
    rtol, atol : float and array_like, optional
        Relative and absolute tolerances. The solver keeps the error estimates
        less than ``atol` + rtol * abs(y)``. Here `rtol` controls a relative
        accuracy (number of correct digits). But if a component of `y` is
        approximately below `atol` then the error only needs to fall within
        the same `atol` threshold, and the number of correct digits is not
        guaranteed. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    method : string, optional
        Integration method to use:

            * 'RK45' (default): Explicit Runge-Kutta method of order 5 with an
              automatic step size control [1]_. A 4-th order accurate quartic
              polynomial is used for the continuous extension [2]_.
            * 'RK23': Explicit Runge-Kutta method of order 3 with an automatic
              step size control [3]_. A 3-th order accurate cubic Hermit
              polynomial is used for the continuous extension.

    events : callable, list of callables or None, optional
        Events to track. Events are defined by functions which take
        a zero value at a point of an event. Each function must have a
        signature ``event(x, y)`` and return float, the solver will find an
        accurate value of ``x`` at which ``event(x, y(x)) = 0`` using a root
        finding algorithm. Additionally each ``event`` function might have
        attributes:

            * terminate : bool, whether to terminate integration if this
              event occurs. Implicitly False if not assigned.
            * direction : float, direction of crossing a zero. If `direction`
              is positive then `event` must go from negative to positive, and
              vice-versa if `direction` is negative. If 0, then either way will
              count. Implicitly 0 if not assigned.

        You can assign attributes like ``event.terminate = True`` to any
        function in Python. If None (default), events won't be tracked.

    Returns
    -------
    Bunch object with the following fields defined:
    sol : PPoly
        Found solution for y as `scipy.interpolate.PPoly` instance, a C1
        continuous spline.
    x : ndarray, shape (n_points,)
        Values of the independent variable at which the solver made steps.
    y : ndarray, shape (n, n_points)
        Solution values at `x`.
    yp : ndarray, shape (n, n_points)
        Solution derivatives at `x`, i.e. ``fun(x, y)``.
    x_events : ndarray, tuple of ndarray or None
        Arrays containing values of x at each corresponding events was
        detected. If `events` contained only 1 event, then `x_events` will
        be ndarray itself. None if `events` was None.
    status : int
        Reason for algorithm termination:

            * 0: The solver failed to reach the interval end or a termination
              event.
            * 1: The solver successfully reached the interval end.
            * 2: A termination event occurred.

    message : string
        Verbal description of the termination reason.
    success : bool
        True if the solver reached the interval end or a termination event
        (``status > 0``).

    References
    ----------
    .. [1] J. R. Dormand, P. J. Prince, "A family of embedded Runge-Kutta
           formulae", Journal of Computational and Applied Mathematics, Vol. 6,
           No. 1, pp. 19-26, 1980.
    .. [2] L. W. Shampine, "Some Practical Runge-Kutta Formulas", Mathematics
           of Computation,, Vol. 46, No. 173, pp. 135-150, 1986.
    .. [3] P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
           Appl. Math. Lett. Vol. 2, No. 4. pp. 321-325, 1989.
    """
    if method not in METHODS:
        raise ValueError("`method` must be one of {}.".format(METHODS))

    ya = np.atleast_1d(ya)
    if ya.ndim != 1:
        raise ValueError("`ya` must be 1-dimensional.")

    a, b = float(x_span[0]), float(x_span[1])
    if a == b:
        raise ValueError("Initial and final `x` must be distinct.")

    def fun_wrapped(x, y):
        return np.asarray(fun(x, y), dtype=float)

    fa = fun_wrapped(a, ya)
    if fa.shape != ya.shape:
        raise ValueError("`fun` return is expected to have shape {}, "
                         "but actually has {}.".format(ya.shape, fa.shape))

    if callable(events):
        events = (events,)

    if events is not None:
        direction = np.empty(len(events))
        is_terminal = np.empty(len(events), dtype=bool)
        for i, event in enumerate(events):
            try:
                is_terminal[i] = event.terminate
            except AttributeError:
                is_terminal[i] = False

            try:
                direction[i] = event.direction
            except AttributeError:
                direction[i] = 0
    else:
        direction = None
        is_terminal = None

    status, sol, xs, ys, fs, x_events = rk(
        fun_wrapped, a, b, ya, fa, rtol, atol, method, events, direction,
        is_terminal)

    return ODEResult(sol=sol, x=xs, y=ys, yp=fs, x_events=x_events,
                     status=status, message=TERMINATION_MESSAGES[status],
                     success=status > 0)
