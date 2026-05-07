"""Oscillatory initial value problem solver based on a Riccati formulation."""

from typing import Any, Callable, Mapping, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .ivp import OdeResult
from . import pyriccaticpp as _ric

"""
Mapping of (float, complex) data types to osc C++ class
"""
_INIT_CLASSES: Mapping[Tuple[bool, bool], Any] = {
    (False, False): _ric.Init_OF64_GF64,
    (True,  False): _ric.Init_OC64_GF64,
    (False, True ): _ric.Init_OF64_GC64,
    (True,  True ): _ric.Init_OC64_GC64,
}


def _validate_time_span(t_span: Sequence[float]) -> Tuple[float, float, float]:
    """Return ``(t0, tf, direction)`` from a 2-element time interval."""
    if len(t_span) != 2:
        raise ValueError("`t_span` must be a 2-element sequence.")
    try:
        t0: float = float(t_span[0])
        tf: float = float(t_span[1])
    except (TypeError, ValueError) as exc:
        raise ValueError("Values in `t_span` must be real numbers.") from exc
    if not (np.isfinite(t0) and np.isfinite(tf)):
        raise ValueError("Values in `t_span` must be finite.")
    direction: float = -1.0 if tf < t0 else 1.0 if tf > t0 else 0.0
    return t0, tf, direction


def _validate_y0(y0: ArrayLike) -> NDArray[np.complex128]:
    """Coerce `y0` to a 1-D complex array of length 2."""
    try:
        y0_arr: NDArray[np.complex128] = np.asarray(y0, dtype=np.complex128)
    except (TypeError, ValueError) as exc:
        raise TypeError("`y0` must be array_like.") from exc
    if y0_arr.ndim != 1:
        raise ValueError("`y0` must be 1-dimensional.")
    if len(y0_arr) != 2:
        raise ValueError("`y0` must have exactly 2 elements: [y, y'].")
    return y0_arr


def _validate_positive_float(name: str, value: Any) -> float:
    """Return ``value`` coerced to a positive float, with `name` in errors."""
    try:
        v: float = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"`{name}` must be a positive float.") from exc
    if not (v > 0):
        raise ValueError(f"`{name}` must be positive.")
    return v


def _validate_epsilon_h(epsilon_h: Optional[float], rtol: float) -> float:
    """Validate ``epsilon_h`` or derive ``max(1e-6, rtol)`` when None."""
    if epsilon_h is None:
        return max(1e-6, rtol)
    return _validate_positive_float("epsilon_h", epsilon_h)


def _validate_init_stepsize(init_stepsize: float, direction: float) -> float:
    """Return ``init_stepsize`` as a non-zero float signed to match ``direction``."""
    try:
        value: float = float(init_stepsize)
    except (TypeError, ValueError) as exc:
        raise TypeError("`init_stepsize` must be a float.") from exc
    if value == 0:
        raise ValueError("`init_stepsize` must be non-zero.")
    if direction != 0.0:
        return float(np.copysign(abs(value), direction))
    return value


def _validate_t_eval(
    t_eval: Optional[ArrayLike],
    t0: float,
    tf: float,
    direction: float,
) -> Optional[NDArray[np.float64]]:
    """Validate `t_eval` and return it as a (possibly empty) array, or None.

    Checks that values are 1-D, lie within ``t_span``, and are strictly
    monotonic in the integration direction.
    """
    if t_eval is None:
        return None
    t_eval_array: NDArray[np.float64] = np.asarray(t_eval, dtype=np.float64)
    if t_eval_array.ndim != 1:
        raise ValueError("`t_eval` must be 1-dimensional.")
    if t_eval_array.size == 0 or direction == 0.0:
        return t_eval_array
    t_eval_min: float = min(t0, tf)
    t_eval_max: float = max(t0, tf)
    if np.any(t_eval_array < t_eval_min) or np.any(t_eval_array > t_eval_max):
        raise ValueError("`t_eval` must be within `t_span`.")
    diffs: NDArray[np.float64] = np.diff(t_eval_array)
    if direction > 0 and np.any(diffs <= 0):
        raise ValueError("`t_eval` must be strictly increasing.")
    if direction < 0 and np.any(diffs >= 0):
        raise ValueError("`t_eval` must be strictly decreasing.")
    return t_eval_array


def _validate_solver_options(
    options: Mapping[str, Any],
) -> Tuple[bool, int, int, int, int]:
    """Validate ``**options`` and return ``(hard_stop, nini, nmax, n, p)``."""
    hard_stop: bool = options.get('hard_stop', False)
    if not isinstance(hard_stop, (bool, np.bool_)):
        raise TypeError("`hard_stop` must be bool.")

    nini: int = options.get('nini', 16)
    nmax: int = options.get('nmax', 32)
    # When only one of n/p is provided, mirror it onto the other.
    n: int = options.get('n', options.get('p', 32))
    p: int = options.get('p', options.get('n', 32))

    int_options: list[int] = []
    for name, value in [('nini', nini), ('nmax', nmax), ('n', n), ('p', p)]:
        try:
            value_int: int = int(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"`{name}` must be an integer.") from exc
        if value_int <= 0:
            raise ValueError(f"`{name}` must be positive.")
        int_options.append(value_int)

    return hard_stop, int_options[0], int_options[1], int_options[2], int_options[3]


def _make_solver_info(
    omega_fun: Callable[[ArrayLike], ArrayLike],
    gamma_fun: Callable[[ArrayLike], ArrayLike],
    t0: float,
    nini: int, nmax: int, n: int, p: int,
) -> Any:
    """Probe omega/gamma at t0 and pick the matching pyriccaticpp Init_* class.

    The riccati core calls both callbacks with scalar floats *and* numpy
    arrays during a single solve. Failing the array call here surfaces a
    clear Python error instead of an opaque pybind11 cast failure later.
    """
    t0_scalar: float = float(t0)
    t0_array: NDArray[np.float64] = np.array([t0_scalar])
    sample_omega = omega_fun(t0_scalar)
    sample_gamma = gamma_fun(t0_scalar)
    for name, fn in (("omega_fun", omega_fun), ("gamma_fun", gamma_fun)):
        try:
            fn(t0_array)
        except Exception as exc:
            raise TypeError(
                f"`{name}` must accept array inputs as well as scalars; "
                f"the riccati core calls it both ways."
            ) from exc
    cls = _INIT_CLASSES[
        (isinstance(sample_omega, complex), isinstance(sample_gamma, complex))
    ]
    return cls(omega_fun, gamma_fun, nini, nmax, n, p)


def solve_ivp_osc(
    omega_fun: Callable[[ArrayLike], ArrayLike],
    gamma_fun: Callable[[ArrayLike], ArrayLike],
    t_span: Sequence[float],
    y0: ArrayLike,
    t_eval: Optional[ArrayLike] = None,
    rtol: float = 1e-3,
    atol: float = 1e-6,
    epsilon_h: Optional[float] = None,
    init_stepsize: float = 0.01,
    **options: Any,
) -> OdeResult:
    """Solve an initial value problem for oscillatory systems.

    This function is a specialized variant of :func:`solve_ivp` for
    oscillatory problems. It integrates second-order oscillatory systems
    using the adaptive Riccati defect correction (ARDC) formulation.

    Parameters
    ----------
    omega_fun, gamma_fun : callable
        Callables describing the oscillation frequency omega(x) and
        friction gamma(x). Both must accept a scalar or array argument
        and return values of the same shape. They can return either
        real or complex values.
    t_span : 2-member sequence of float
        Interval of integration ``(t0, tf)``. Both entries must be
        finite, real numbers.
    y0 : array_like, shape (2,)
        Initial conditions as ``[y(t0), y'(t0)]``. Must be a 2-element
        array containing the initial value and its derivative. Can be
        complex-valued.
    t_eval : array_like or None, optional
        Times at which to store the computed solution. When None (the
        default), no dense interpolation is performed and the returned
        ``t``/``y``/``ydot`` are the adaptive step endpoints chosen by
        the solver. When provided, the values must be strictly monotonic
        in the integration direction and lie within `t_span`; the
        riccati core then interpolates the solution at these points
        using its per-step Chebyshev coefficients in a single pass, and
        the adaptive step grid is *not* surfaced -- ``t``/``y``/``ydot``
        contain only the values at `t_eval`. The number of adaptive
        steps the solver actually took is still recoverable as
        ``len(result.extra['successes'])``.
    rtol, atol : float, optional
        Relative and absolute tolerances. The Riccati core takes a single
        relative tolerance `eps`; the two are folded together via
        ``eps = max(rtol + atol / max(|y0[0]|, |y0[1]|, 1.0), 1e-15)``,
        approximating the standard ``|err| < rtol*|y| + atol`` mixed
        tolerance evaluated at ``t0``.
    epsilon_h : float, optional
        Tolerance for stepsize selection. Controls how accurately the
        frequency and friction functions are interpolated. If not provided,
        defaults to `max(1e-6, rtol)`.
    init_stepsize : float, optional
        Initial stepsize for integration. Default is 0.01. The
        magnitude is what the user controls; the sign is taken from
        the integration direction (``tf - t0``), so passing a positive
        value also works for backward integration.
    **options
        Advanced solver options:

        - nini : int, default 16
            Minimum number of Chebyshev nodes
        - nmax : int, default 32
            Maximum number of Chebyshev nodes
        - n : int, default 32
            Number of Chebyshev nodes for collocation steps
        - p : int, default same as n (32)
            Number of Chebyshev nodes for Riccati steps
        - hard_stop : bool, default False
            If True, force the solver to land exactly on `tf` (potentially
            with a smaller last step).

    Returns
    -------
    OdeResult
        Object with fields:

        - t : ndarray, shape (n_points,), dtype float
            Time points. Equals `t_eval` when provided; otherwise the
            adaptive step endpoints chosen by the solver. The two cases
            return *different* arrays from different sources -- when
            `t_eval` is given the adaptive step grid is discarded and
            only the user-requested points are surfaced.
        - y : ndarray, shape (n_points,), dtype complex
            Solution values at `t`. When `t_eval` is given these are
            obtained via the riccati core's per-step Chebyshev
            interpolation; otherwise they are the endpoint values
            recorded at each adaptive step.
        - ydot : ndarray, shape (n_points,), dtype complex
            Derivative values at `t`, computed the same way as ``y``.
        - status : int
            Reason for algorithm termination (0 = success).
        - message : str
            Verbal description of status.
        - success : bool
            True if solver reached end of interval.
        - nfev : int
            Number of omega_fun/gamma_fun evaluations. The Riccati core
            does not currently track callback counts, so this is always 0.
        - extra : dict
            Riccati diagnostics: `successes`, `phases`, `steptypes`.

    Raises
    ------
    TypeError
        If inputs have incorrect types (e.g., non-callable functions, non-float
        tolerances, non-boolean flags).
    ValueError
        If `t_span` or `t_eval` are invalid, tolerances are non-positive, or
        `y0` has the wrong shape.

    .. versionadded:: 1.18.0

    Notes
    -----
    The riccati solver is specialized for second-order ODEs of the form
    ``y'' + 2*gamma(x)*y' + omega(x)**2 * y = 0``. It uses adaptive
    Chebyshev spectral collocation and can efficiently handle both
    oscillatory and non-oscillatory regions.

    Examples
    --------
    Solve the Airy equation ``y'' + x*y = 0`` on ``x > 0``, whose
    solutions are ``Ai(-x)`` and ``Bi(-x)``. In the Riccati form
    ``y'' + 2*gamma*y' + omega**2 * y = 0``, this corresponds to
    ``omega(x) = sqrt(x)`` and ``gamma(x) = 0``:

    >>> from scipy.integrate import solve_ivp_osc
    >>> from scipy.special import airy
    >>> import numpy as np
    >>> omega_fun = lambda x: np.sqrt(x)
    >>> gamma_fun = lambda x: np.zeros_like(x)
    >>> x0, xf = 1.0, 100.0
    >>> Ai0, Aip0, Bi0, Bip0 = airy(-x0)
    >>> y0 = [Ai0 + 1j*Bi0, -Aip0 - 1j*Bip0]  # y(x0), y'(x0)
    >>> result = solve_ivp_osc(omega_fun, gamma_fun, (x0, xf), y0)
    """
    if not callable(omega_fun):
        raise TypeError("`omega_fun` must be callable.")
    if not callable(gamma_fun):
        raise TypeError("`gamma_fun` must be callable.")
    t0: float
    tf: float
    direction: float
    t0, tf, direction = _validate_time_span(t_span)
    y0_array: NDArray[np.complex128] = _validate_y0(y0)
    rtol_: float = _validate_positive_float("rtol", rtol)
    atol_: float = _validate_positive_float("atol", atol)
    # Fold (rtol, atol) into a single relative tolerance using the y0 scale,
    # approximating the standard mixed tolerance |err| < rtol*|y| + atol.
    y_scale: float = max(abs(y0_array[0]), abs(y0_array[1]), 1.0)
    eps: float = max(rtol_ + atol_ / y_scale, 1e-15)
    epsilon_h_value: float = _validate_epsilon_h(epsilon_h, rtol_)
    init_stepsize_value: float = _validate_init_stepsize(init_stepsize, direction)
    t_eval_array: Optional[NDArray[np.float64]] = _validate_t_eval(
        t_eval, t0, tf, direction,
    )
    hard_stop: bool
    nini: int
    nmax: int
    n: int
    p: int
    hard_stop, nini, nmax, n, p = _validate_solver_options(options)

    info = _make_solver_info(omega_fun, gamma_fun, t0, nini, nmax, n, p)
    x_eval_arg = (t_eval_array
                  if t_eval_array is not None and t_eval_array.size > 0
                  else None)
    try:
        out = _ric.evolve(
            info, t0, tf,
            complex(y0_array[0]), complex(y0_array[1]),
            eps, epsilon_h_value, init_stepsize_value,
            x_eval_arg, hard_stop, _ric.LogLevel.ERROR,
        )
    except Exception as e:
        fail = OdeResult(
            t=np.array([t0]),
            y=y0_array[0:1],
            nfev=0,
            status=-1,
            message=f"Riccati solver failed: {e}",
            success=False,
        )
        fail.ydot = y0_array[1:2]
        fail.extra = {
            "successes": np.empty(0, dtype=int),
            "phases": np.empty(0, dtype=float),
            "steptypes": np.empty(0, dtype=int),
        }
        return fail
    xs, ys, dys, success_out, phase_out, steptype_out, yeval, dyeval, _ = out

    status: int
    message: str
    success: bool
    if success_out and success_out[-1] == 1:
        status = 0
        message = "The solver successfully reached the end of the integration interval."
        success = True
    else:
        status = -1
        message = "The solver did not reach the end of the integration interval."
        success = False

    if t_eval_array is not None:
        t_out: NDArray[np.float64] = t_eval_array
        y_out: NDArray[np.complex128] = np.asarray(yeval, dtype=np.complex128)
        ydot_out: NDArray[np.complex128] = np.asarray(dyeval, dtype=np.complex128)
    else:
        t_out = np.asarray(xs, dtype=np.float64)
        y_out = np.asarray(ys, dtype=np.complex128)
        ydot_out = np.asarray(dys, dtype=np.complex128)

    result: OdeResult = OdeResult(
        t=t_out,
        y=y_out,
        nfev=0,
        status=status,
        message=message,
        success=success,
    )
    result.ydot = ydot_out
    result.extra = {
        "successes": np.asarray(success_out),
        "phases": np.asarray(phase_out),
        "steptypes": np.asarray(steptype_out),
    }
    return result
