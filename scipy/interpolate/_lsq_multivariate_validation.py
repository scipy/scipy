"""Input validation helpers for least-squares multivariate splines."""

from __future__ import annotations

from dataclasses import dataclass
from operator import index

import numpy as np


@dataclass
class _SplineInput:
    """Validated and normalized constructor inputs."""

    x: object
    y: object
    t: object
    w: object
    bbox: object
    k: object
    eps: object
    check_finite: bool
    sparse: bool
    smoothing: float
    penalty_order: int
    unsupported: str


def _validate_input(
    x,
    y,
    t,
    w,
    bbox,
    k,
    eps,
    check_finite,
    sparse,
    smoothing,
    penalty_order,
    unsupported,
):
    """Validate and normalize constructor inputs."""
    y = np.asarray(y, dtype=float)

    if check_finite and not np.all(np.isfinite(y)):
        raise ValueError("y must contain only finite values.")

    if y.ndim not in (1, 2):
        raise ValueError("y must be a 1D or 2D array.")

    if y.ndim == 2 and y.shape[1] == 0:
        raise ValueError("y must contain at least one output dimension.")

    n_samples = y.shape[0]
    if n_samples == 0:
        raise ValueError("x and y must contain at least one sample.")

    x = np.asarray(x, dtype=float)

    if check_finite and not np.all(np.isfinite(x)):
        raise ValueError("x must contain only finite values.")

    if x.ndim == 1:
        x = x[:, None]
    elif x.ndim == 2:
        if x.shape[0] == n_samples:
            pass
        elif x.shape[1] == n_samples:
            x = x.T
        else:
            raise ValueError(
                "x must have shape (n_samples, n_dimensions) or "
                "(n_dimensions, n_samples)."
            )
    else:
        raise ValueError("x must be a 1D or 2D array.")

    if y.shape[0] != x.shape[0]:
        raise ValueError("x and y must contain the same number of samples.")

    n_samples, ndim = x.shape
    if ndim == 0:
        raise ValueError("x must contain at least one input dimension.")

    w = _validate_weights(w, n_samples, check_finite=check_finite)
    degrees = _normalize_degrees(k, ndim)
    bbox = _normalize_bbox(bbox, x, ndim, check_finite=check_finite)
    knot_vectors = _normalize_knot_input(t, x=x, bbox=bbox, ndim=ndim)
    _validate_knot_vectors(knot_vectors, bbox, check_finite=check_finite)
    eps = _validate_eps(eps)
    smoothing = _validate_smoothing(smoothing)
    penalty_order = _validate_penalty_order(penalty_order)
    unsupported = _validate_unsupported(unsupported)

    if not isinstance(sparse, (bool, np.bool_)):
        raise ValueError("sparse must be True or False.")

    return _SplineInput(
        x=x,
        y=y,
        t=knot_vectors,
        w=w,
        bbox=bbox,
        k=degrees,
        eps=eps,
        check_finite=check_finite,
        sparse=bool(sparse),
        smoothing=smoothing,
        penalty_order=penalty_order,
        unsupported=unsupported,
    )


def _validate_weights(w, n_samples, check_finite):
    """Validate weights and fill the default unweighted case."""
    if w is None:
        return np.ones(n_samples, dtype=float)

    w = np.asarray(w, dtype=float)

    if check_finite and not np.all(np.isfinite(w)):
        raise ValueError("w must contain only finite values.")

    if w.ndim != 1:
        raise ValueError("w must be a 1D array.")

    if w.shape[0] != n_samples:
        raise ValueError("w must have the same length as y.")

    if np.any(w <= 0):
        raise ValueError("w must contain only positive values.")

    return w


def _normalize_degrees(k, ndim):
    """Normalize spline degree input to one integer per dimension."""
    if np.isscalar(k):
        try:
            degree = index(k)
        except TypeError as exc:
            raise ValueError(
                "k must be an integer or a sequence of integers."
            ) from exc

        degrees = (degree,) * ndim
    else:
        try:
            degrees = tuple(index(ki) for ki in k)
        except TypeError as exc:
            raise ValueError(
                "k must be an integer or a sequence of integers."
            ) from exc

    if len(degrees) != ndim:
        raise ValueError("k must be a scalar or have one value per dimension.")

    if any(ki < 1 for ki in degrees):
        raise ValueError("spline degrees must be positive integers.")

    return degrees


def _normalize_bbox(bbox, x, ndim, check_finite):
    """Normalize bbox to shape (n_dimensions, 2)."""
    if bbox is None:
        bbox = np.column_stack((np.min(x, axis=0), np.max(x, axis=0)))
    else:
        bbox = np.asarray(bbox, dtype=float)

        if bbox.shape == (2,) and ndim == 1:
            bbox = bbox[None, :]
        elif bbox.shape == (2 * ndim,):
            bbox = bbox.reshape(ndim, 2)

        if bbox.shape != (ndim, 2):
            raise ValueError("bbox must have shape (n_dimensions, 2).")

    if check_finite and not np.all(np.isfinite(bbox)):
        raise ValueError("bbox must contain only finite values.")

    if np.any(bbox[:, 0] >= bbox[:, 1]):
        raise ValueError("each bbox lower bound must be less than its upper bound.")

    if np.any(x < bbox[:, 0]) or np.any(x > bbox[:, 1]):
        raise ValueError("all x values must lie inside bbox.")

    return bbox


def _normalize_knot_input(t, x, bbox, ndim):
    """Normalize explicit or automatic knot input."""
    if _is_scalar_like(t):
        raw_knots = (t,) * ndim
    elif ndim == 1 and not _is_sequence_of_sequences(t):
        raw_knots = (t,)
    else:
        try:
            raw_knots = tuple(t)
        except TypeError as exc:
            raise ValueError(
                "t must be an integer or contain one knot specification "
                "per dimension."
            ) from exc

        if len(raw_knots) != ndim:
            raise ValueError("t must contain one knot vector per dimension.")

    knot_vectors = []
    for axis, knot_spec in enumerate(raw_knots):
        if _is_scalar_like(knot_spec):
            try:
                n_pieces = index(knot_spec)
            except TypeError as exc:
                raise ValueError(
                    "automatic knot counts in t must be integers."
                ) from exc

            knots = _automatic_interior_knots(x[:, axis], bbox[axis], n_pieces)
        else:
            knots = np.asarray(knot_spec, dtype=float)

        knot_vectors.append(knots)

    return tuple(knot_vectors)


def _automatic_interior_knots(x, bbox, n_pieces):
    """Choose interior knots from an intended number of polynomial pieces."""
    if n_pieces < 1:
        raise ValueError("automatic knot counts in t must be positive.")

    n_interior = n_pieces - 1
    if n_interior == 0:
        return np.empty(0, dtype=float)

    lower, upper = bbox
    probabilities = np.linspace(0.0, 1.0, n_pieces + 1)[1:-1]
    knots = np.quantile(x, probabilities)
    knots = knots[(lower < knots) & (knots < upper)]
    knots = np.unique(knots)

    if knots.size != n_interior:
        knots = np.linspace(lower, upper, n_pieces + 1)[1:-1]

    return np.asarray(knots, dtype=float)


def _validate_knot_vectors(knot_vectors, bbox, check_finite):
    """Validate normalized interior knot vectors."""
    for axis, knots in enumerate(knot_vectors):
        if knots.ndim != 1:
            raise ValueError("each knot vector in t must be 1D.")

        if check_finite and not np.all(np.isfinite(knots)):
            raise ValueError("t must contain only finite knot values.")

        if knots.size and np.any(np.diff(knots) <= 0):
            raise ValueError("each knot vector in t must be strictly increasing.")

        lower, upper = bbox[axis]
        if knots.size and (np.any(knots <= lower) or np.any(knots >= upper)):
            raise ValueError("interior knots must lie strictly inside bbox.")


def _validate_eps(eps):
    """Validate least-squares rank threshold."""
    if eps is None:
        return eps

    if not np.isscalar(eps):
        raise ValueError("eps must be a scalar.")

    eps = float(eps)
    if not np.isfinite(eps):
        raise ValueError("eps must be finite.")

    if not 0.0 < eps < 1.0:
        raise ValueError("eps must lie in the open interval (0, 1).")

    return eps


def _validate_smoothing(smoothing):
    """Validate roughness-penalty smoothing strength."""
    if not np.isscalar(smoothing):
        raise ValueError("smoothing must be a scalar.")

    smoothing = float(smoothing)
    if not np.isfinite(smoothing):
        raise ValueError("smoothing must be finite.")

    if smoothing < 0.0:
        raise ValueError("smoothing must be nonnegative.")

    return smoothing


def _validate_penalty_order(penalty_order):
    """Validate coefficient finite-difference penalty order."""
    try:
        penalty_order = index(penalty_order)
    except TypeError as exc:
        raise ValueError("penalty_order must be an integer.") from exc

    if penalty_order < 1:
        raise ValueError("penalty_order must be a positive integer.")

    return penalty_order


def _validate_unsupported(unsupported):
    """Validate unsupported basis-function handling."""
    if unsupported not in {"raise", "warn", "ignore"}:
        raise ValueError("unsupported must be 'raise', 'warn', or 'ignore'.")

    return unsupported


def _validate_eval_input(x, ndim):
    """Validate and normalize evaluation coordinates."""
    was_scalar = np.ndim(x) == 0
    x = np.asarray(x, dtype=float)

    if x.ndim == 0:
        if ndim != 1:
            raise ValueError("scalar evaluation is only valid for 1D splines.")
        x = x[None, None]
    elif x.ndim == 1:
        if ndim == 1:
            x = x[:, None]
        elif x.shape[0] == ndim:
            x = x[None, :]
        else:
            raise ValueError("x must have shape (n_eval, n_dimensions).")
    elif x.ndim == 2:
        if x.shape[1] != ndim:
            raise ValueError("x must have shape (n_eval, n_dimensions).")
    else:
        raise ValueError("x must be a scalar, 1D array, or 2D array.")

    return x, was_scalar


def _validate_eval_domain(x, bbox):
    """Validate that evaluation coordinates are inside the spline domain."""
    if np.any(x < bbox[:, 0]) or np.any(x > bbox[:, 1]):
        raise ValueError("evaluation coordinates must lie inside bbox.")


def _validate_derivative_orders(nu, ndim):
    """Validate derivative orders for evaluation."""
    if nu is None:
        return (0,) * ndim

    if np.isscalar(nu):
        try:
            order = index(nu)
        except TypeError as exc:
            raise ValueError("nu must contain integer derivative orders.") from exc

        orders = (order,) * ndim
    else:
        try:
            orders = tuple(index(order) for order in nu)
        except TypeError as exc:
            raise ValueError("nu must contain integer derivative orders.") from exc

        if len(orders) != ndim:
            raise ValueError("nu must be a scalar or have one value per dimension.")

    if any(order < 0 for order in orders):
        raise ValueError("derivative orders in nu must be nonnegative.")

    return orders


def _is_sequence_of_sequences(value):
    """Return True when value looks like a sequence of knot vectors."""
    if isinstance(value, (str, bytes)):
        return False

    try:
        items = list(value)
    except TypeError:
        return False

    if not items:
        return False

    return all(np.ndim(item) > 0 for item in items)


def _is_scalar_like(value):
    """Return True for scalar-like values without coercing ragged sequences."""
    try:
        return np.ndim(value) == 0
    except ValueError:
        return False
