"""Grid-data conversion helpers for multivariate splines."""

from __future__ import annotations

import numpy as np


def _grid_to_scattered(axes, y, w=None):
    """Convert tensor-grid coordinates and values to scattered rows."""
    try:
        axes = tuple(np.asarray(axis, dtype=float) for axis in axes)
    except TypeError as exc:
        raise ValueError("axes must contain one coordinate vector per axis.") from exc

    if not axes:
        raise ValueError("axes must contain at least one coordinate vector.")

    for axis in axes:
        if axis.ndim != 1:
            raise ValueError("each grid axis must be a 1D array.")
        if axis.size == 0:
            raise ValueError("grid axes must be nonempty.")

    grid_shape = tuple(axis.size for axis in axes)
    mesh = np.meshgrid(*axes, indexing="ij")
    x = np.column_stack([axis_grid.ravel() for axis_grid in mesh])

    y = _reshape_grid_values(y, grid_shape)
    if w is not None:
        w = _reshape_grid_weights(w, grid_shape, n_samples=x.shape[0])

    return x, y, w


def _reshape_grid_values(y, grid_shape):
    """Reshape scalar or vector-valued grid data to sample rows."""
    y = np.asarray(y, dtype=float)
    ndim = len(grid_shape)

    if y.shape == grid_shape:
        return y.ravel()

    if y.ndim == ndim + 1 and y.shape[:-1] == grid_shape:
        return y.reshape(-1, y.shape[-1])

    if y.ndim == ndim + 1 and y.shape[1:] == grid_shape:
        return np.moveaxis(y, 0, -1).reshape(-1, y.shape[0])

    raise ValueError(
        "grid y must have shape grid_shape, grid_shape + (n_outputs,), "
        "or (n_outputs,) + grid_shape."
    )


def _reshape_grid_weights(w, grid_shape, n_samples):
    """Reshape grid weights to sample rows."""
    w = np.asarray(w, dtype=float)
    if w.shape == grid_shape:
        return w.ravel()

    if w.shape != (n_samples,):
        raise ValueError("grid weights must match the grid shape.")

    return w
