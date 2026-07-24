"""Roughness penalty matrices for tensor-product spline coefficients."""

from __future__ import annotations

from math import comb

import numpy as np


def _roughness_penalty_matrix(coefficient_shape, penalty_order):
    """Build a sparse finite-difference roughness penalty matrix."""
    sparse = _require_scipy_sparse()
    coefficient_shape = tuple(coefficient_shape)
    n_coeffs = int(np.prod(coefficient_shape))
    axis_penalties = []

    for axis, n_axis in enumerate(coefficient_shape):
        if n_axis <= penalty_order:
            continue

        axis_penalties.append(
            _axis_difference_penalty(
                coefficient_shape,
                axis,
                penalty_order,
                sparse,
            )
        )

    if not axis_penalties:
        return sparse.csr_matrix((0, n_coeffs))

    return sparse.vstack(axis_penalties, format="csr")


def _axis_difference_penalty(coefficient_shape, axis, penalty_order, sparse):
    """Build the finite-difference penalty for one coefficient axis."""
    operators = [
        sparse.identity(size, format="csr")
        for size in coefficient_shape
    ]
    operators[axis] = _difference_matrix(
        coefficient_shape[axis],
        penalty_order,
        sparse,
    )

    penalty = operators[0]
    for operator in operators[1:]:
        penalty = sparse.kron(penalty, operator, format="csr")

    return penalty


def _difference_matrix(n, penalty_order, sparse):
    """Build a 1D finite-difference matrix."""
    coefficients = [
        (-1) ** (penalty_order - offset) * comb(penalty_order, offset)
        for offset in range(penalty_order + 1)
    ]
    offsets = np.arange(penalty_order + 1)

    return sparse.diags(
        coefficients,
        offsets,
        shape=(n - penalty_order, n),
        format="csr",
    )


def _require_scipy_sparse():
    """Import scipy.sparse for smoothing penalty construction."""
    try:
        from scipy import sparse
    except ImportError as exc:
        raise ImportError("smoothing > 0 requires scipy to be installed.") from exc

    return sparse
