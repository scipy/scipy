"""Least-squares solvers for spline coefficients."""

from __future__ import annotations

import warnings

import numpy as np

from ._lsq_multivariate_penalty import _roughness_penalty_matrix


def _solve_lsq(
    design_matrix,
    y,
    w=None,
    eps=None,
    smoothing=0.0,
    coefficient_shape=None,
    penalty_order=2,
    allow_rank_deficient=False,
    warn_rank_deficient=False,
):
    """Solve the weighted least-squares system for spline coefficients."""
    is_sparse = _is_sparse_matrix(design_matrix)

    if w is None:
        weighted_matrix = design_matrix
        weighted_y = y
    else:
        weighted_matrix = _apply_weights_to_matrix(
            design_matrix,
            w,
            is_sparse=is_sparse,
        )
        weighted_y = _apply_weights_to_y(y, w)

    if smoothing > 0.0:
        weighted_matrix, weighted_y = _augment_with_penalty(
            weighted_matrix,
            weighted_y,
            smoothing=smoothing,
            coefficient_shape=coefficient_shape,
            penalty_order=penalty_order,
        )
        is_sparse = _is_sparse_matrix(weighted_matrix)

    if is_sparse:
        coeffs = _solve_sparse_lsq(weighted_matrix, weighted_y, eps=eps)
    else:
        coeffs, _, rank, _ = np.linalg.lstsq(
            weighted_matrix,
            weighted_y,
            rcond=eps,
        )
        if (
            smoothing == 0.0
            and rank < weighted_matrix.shape[1]
            and not allow_rank_deficient
        ):
            raise ValueError(
                "least-squares system is rank deficient; adjust knots, "
                "sample locations, or use smoothing > 0."
            )
        if (
            smoothing == 0.0
            and rank < weighted_matrix.shape[1]
            and allow_rank_deficient
            and warn_rank_deficient
        ):
            warnings.warn(
                "least-squares system is rank deficient; continuing with "
                "the minimal-norm least-squares solution.",
                RuntimeWarning,
                stacklevel=2,
            )

    residual = _weighted_residual_sum(design_matrix, coeffs, y, w)

    return coeffs, residual


def _augment_with_penalty(
    weighted_matrix,
    weighted_y,
    smoothing,
    coefficient_shape,
    penalty_order,
):
    """Append roughness penalty rows to a weighted least-squares system."""
    if coefficient_shape is None:
        raise ValueError("coefficient_shape is required when smoothing > 0.")

    sparse = _require_scipy_sparse()
    penalty = _roughness_penalty_matrix(coefficient_shape, penalty_order)
    penalty = np.sqrt(smoothing) * penalty

    if penalty.shape[0] == 0:
        return weighted_matrix, weighted_y

    weighted_matrix = sparse.vstack(
        (
            sparse.csr_matrix(weighted_matrix),
            penalty,
        ),
        format="csr",
    )

    if weighted_y.ndim == 1:
        zero_targets = np.zeros(penalty.shape[0], dtype=weighted_y.dtype)
        weighted_y = np.concatenate((weighted_y, zero_targets))
    else:
        zero_targets = np.zeros(
            (penalty.shape[0], weighted_y.shape[1]),
            dtype=weighted_y.dtype,
        )
        weighted_y = np.vstack((weighted_y, zero_targets))

    return weighted_matrix, weighted_y


def _apply_weights_to_matrix(design_matrix, w, is_sparse):
    """Apply row weights to a dense or sparse design matrix."""
    if is_sparse:
        return design_matrix.multiply(w[:, None])

    return design_matrix * w[:, None]


def _apply_weights_to_y(y, w):
    """Apply sample weights to scalar or vector-valued observations."""
    if y.ndim == 1:
        return y * w

    return y * w[:, None]


def _weighted_residual_sum(design_matrix, coeffs, y, w):
    """Compute the weighted sum of squared residuals."""
    residual_values = design_matrix @ coeffs - y
    if w is None:
        weighted_residuals = residual_values
    elif y.ndim == 1:
        weighted_residuals = w * residual_values
    else:
        weighted_residuals = w[:, None] * residual_values

    return np.sum(weighted_residuals**2)


def _solve_sparse_lsq(weighted_matrix, weighted_y, eps=None):
    """Solve a sparse least-squares problem column by column."""
    sparse_linalg = _require_scipy_sparse_linalg()
    tolerance = 1e-12 if eps is None else eps

    if weighted_y.ndim == 1:
        result = sparse_linalg.lsqr(
            weighted_matrix,
            weighted_y,
            atol=tolerance,
            btol=tolerance,
        )
        return result[0]

    coeffs = [
        sparse_linalg.lsqr(
            weighted_matrix,
            weighted_y[:, output],
            atol=tolerance,
            btol=tolerance,
        )[0]
        for output in range(weighted_y.shape[1])
    ]
    return np.column_stack(coeffs)


def _is_sparse_matrix(matrix):
    """Return True if matrix is a SciPy sparse matrix."""
    try:
        from scipy import sparse
    except ImportError:
        return False

    return sparse.issparse(matrix)


def _require_scipy_sparse_linalg():
    """Import scipy.sparse.linalg only when sparse solving is requested."""
    try:
        from scipy.sparse import linalg
    except ImportError as exc:
        raise ImportError("sparse=True requires scipy to be installed.") from exc

    return linalg


def _require_scipy_sparse():
    """Import scipy.sparse for penalized least-squares systems."""
    try:
        from scipy import sparse
    except ImportError as exc:
        raise ImportError("smoothing > 0 requires scipy to be installed.") from exc

    return sparse
