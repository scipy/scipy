"""Generic interface for initial value problem solvers."""
from __future__ import division, print_function, absolute_import

from warnings import warn
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from .common import EPS, ODEResult
from .rk import rk
from .radau import radau


METHODS = ['RK23', 'RK45', 'Radau']


TERMINATION_MESSAGES = {
    -1: "Required step size became too small.",
    0: "The solver successfully reached the interval end.",
    1: "A termination event occurred."
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


def prepare_events(events):
    if callable(events):
        events = (events,)

    if events is not None:
        is_terminal = np.empty(len(events), dtype=bool)
        direction = np.empty(len(events))
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
        is_terminal = None
        direction = None

    return events, is_terminal, direction


def solve_ivp(fun, x_span, ya, rtol=1e-3, atol=1e-6, method='RK45', M=None,
              max_step=None, jac=None, events=None):
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
            * 'Radau': Implicit Runge-Kutta method of Radau IIA family of
              order 5 [4]_. A 5-th order accurate cubic polynomial is available
              naturally as the method can be viewed as a collocation method.

        You should use 'RK45' or 'RK23' methods for non-stiff problems and
        'Radau' for stiff problems [5]_. If not sure, first try to run 'RK45'
        and if it does unusual many iterations or diverges then your problem
        is likely to be stiff and you should use 'Radau'.
    max_step : float or None, optional
        Maximum allowed step size. If None, a step size is selected to be 0.1
        of the length of `x_span`.
    jac : array_like, callable or None, optional
        Jacobian matrix of the right-hand side of the system with respect to
        `y`, required only by 'Radau' method. The Jacobian matrix has shape
        (n, n) and its element (i, j) is equal to ``d f_i / d y_j``.
        There are 3 ways to define the Jacobian:

            * If array_like, then the Jacobian is assumed to be constant.
            * If callable, then the Jacobian is assumed to depend on both
              x and y, and will be called as ``jac(x, y)`` as necessary.
            * If None (default), then the Jacobian will be approximated by
              finite differences.

        It is generally recommended to provided the Jacobian rather then
        relying on finite difference approximation.
    events : callable, list of callables or None, optional
        Events to track. Events are defined by functions which take
        a zero value at a point of an event. Each function must have a
        signature ``event(x, y)`` and return float, the solver will find an
        accurate value of ``x`` at which ``event(x, y(x)) = 0`` using a root
        finding algorithm. Additionally each ``event`` function might have
        attributes:

            * terminate: bool, whether to terminate integration if this
              event occurs. Implicitly False if not assigned.
            * direction: float, direction of crossing a zero. If `direction`
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

            * -1: Required step size became too small.
            * 0: The solver successfully reached the interval end.
            * 1: A termination event occurred.

    message : string
        Verbal description of the termination reason.
    success : bool
        True if the solver reached the interval end or a termination event
        (``status >= 0``).

    References
    ----------
    .. [1] J. R. Dormand, P. J. Prince, "A family of embedded Runge-Kutta
           formulae", Journal of Computational and Applied Mathematics, Vol. 6,
           No. 1, pp. 19-26, 1980.
    .. [2] L. W. Shampine, "Some Practical Runge-Kutta Formulas", Mathematics
           of Computation,, Vol. 46, No. 173, pp. 135-150, 1986.
    .. [3] P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
           Appl. Math. Lett. Vol. 2, No. 4. pp. 321-325, 1989.
    .. [4] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II:
           Stiff and Differential-Algebraic Problems", Sec. IV.8.
    .. [5] `Stiff equation <https://en.wikipedia.org/wiki/Stiff_equation>`_ on
           Wikipedia.
    """
    if method not in METHODS:
        raise ValueError("`method` must be one of {}.".format(METHODS))

    a, b = float(x_span[0]), float(x_span[1])
    if a == b:
        raise ValueError("Initial and final `x` must be distinct.")

    ya = np.atleast_1d(ya)
    if ya.ndim != 1:
        raise ValueError("`ya` must be 1-dimensional.")
    n = ya.shape[0]

    def fun_wrapped(x, y):
        return np.asarray(fun(x, y), dtype=float)

    fa = fun_wrapped(a, ya)
    if fa.shape != ya.shape:
        raise ValueError("`fun` return is expected to have shape {}, "
                         "but actually has {}.".format(ya.shape, fa.shape))

    if method in ['RK23', 'RK45']:
        if M is not None:
            raise ValueError("`M` is supported only by 'Radau' method.")
    elif method == 'Radau':
        if jac is None:
            def jac_wrapped(x, y):
                return approx_derivative(lambda z: fun(x, z),
                                         y, method='2-point')
            Ja = jac_wrapped(a, ya)
        elif callable(jac):
            def jac_wrapped(x, y):
                return np.asarray(jac(x, y), dtype=float)
            Ja = jac_wrapped(a, ya)
            if Ja.shape != (n, n):
                raise ValueError("`jac` return is expected to have shape {}, "
                                 "but actually has {}."
                                 .format((n, n), Ja.shape))
        else:
            Ja = np.asarray(jac, dtype=float)
            if Ja.shape != (n, n):
                raise ValueError("`jac` is expected to have shape {}, but "
                                 "actually has {}.".format((n, n), Ja.shape))
            jac_wrapped = None

    events, is_terminal, direction = prepare_events(events)

    if max_step is None:
        max_step = 0.1 * np.abs(b - a)

    if method in ['RK23', 'RK45']:
        status, sol, xs, ys, fs, x_events = rk(
            fun_wrapped, a, b, ya, fa, rtol, atol, max_step, method,
            events, is_terminal, direction)
    elif method == 'Radau':
        status, sol, xs, ys, fs, x_events = radau(
            fun, jac_wrapped, a, b, ya, fa, Ja, rtol, atol, max_step,
            events, is_terminal, direction)

    return ODEResult(sol=sol, x=xs, y=ys, yp=fs, x_events=x_events,
                     status=status, message=TERMINATION_MESSAGES[status],
                     success=status >= 0)
