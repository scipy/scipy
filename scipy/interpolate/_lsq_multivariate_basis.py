"""Tensor-product B-spline basis construction."""

from __future__ import annotations

import itertools

import numpy as np


def _construct_full_knots(t, bbox, degrees):
    """Construct full knot vectors from interior knots and boundaries."""
    full_knots = []

    for axis, interior_knots in enumerate(t):
        degree = degrees[axis]
        lower, upper = bbox[axis]

        knots = np.concatenate(
            (
                np.repeat(lower, degree + 1),
                interior_knots,
                np.repeat(upper, degree + 1),
            )
        )
        full_knots.append(knots)

    return tuple(full_knots)


def _basis_counts_from_knots(knots, degrees):
    """Return the number of basis functions along each input axis."""
    return tuple(
        len(axis_knots) - degree - 1
        for axis_knots, degree in zip(knots, degrees)
    )


def _bspline_basis_axis(x, knots, degree, basis_index, derivative_order=0):
    """Evaluate one 1D B-spline basis function or derivative."""
    matrix = _basis_matrix_axis(
        x,
        knots,
        degree,
        derivative_order=derivative_order,
    )
    return matrix[:, basis_index]


def _design_matrix_axis(x, knots, degree, derivative_order=0):
    """Build the 1D basis matrix for one input axis."""
    return _basis_matrix_axis(
        x,
        knots,
        degree,
        derivative_order=derivative_order,
    )


def _basis_matrix_axis(x, knots, degree, derivative_order=0):
    """Evaluate all basis functions for one axis."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    knots = np.asarray(knots, dtype=float)
    n_basis = len(knots) - degree - 1

    if derivative_order > degree:
        return np.zeros((x.size, n_basis), dtype=float)

    if derivative_order > 0:
        return _basis_derivative_matrix_axis(
            x,
            knots,
            degree,
            derivative_order,
        )

    return _basis_value_matrix_axis(x, knots, degree)


def _basis_value_matrix_axis(x, knots, degree):
    """Evaluate all non-derivative basis functions for one axis."""
    values = _degree_zero_basis_matrix(x, knots)

    for current_degree in range(1, degree + 1):
        values = _raise_basis_degree(x, knots, values, current_degree)

    endpoint = x == knots[-1]
    if np.any(endpoint):
        values[endpoint, :] = 0.0
        values[endpoint, -1] = 1.0

    return values


def _degree_zero_basis_matrix(x, knots):
    """Evaluate all degree-zero basis functions."""
    n_basis = len(knots) - 1
    values = np.zeros((x.size, n_basis), dtype=float)

    for basis_index in range(n_basis):
        left = knots[basis_index]
        right = knots[basis_index + 1]
        values[:, basis_index] = (left <= x) & (x < right)

    endpoint = x == knots[-1]
    if np.any(endpoint):
        values[endpoint, -1] = 1.0

    return values


def _raise_basis_degree(x, knots, lower_values, degree):
    """Raise a basis matrix by one degree with Cox-de Boor recursion."""
    n_basis = len(knots) - degree - 1
    values = np.zeros((x.size, n_basis), dtype=float)

    for basis_index in range(n_basis):
        left_denominator = knots[basis_index + degree] - knots[basis_index]
        if left_denominator != 0:
            values[:, basis_index] += (
                (x - knots[basis_index])
                / left_denominator
                * lower_values[:, basis_index]
            )

        right_denominator = (
            knots[basis_index + degree + 1] - knots[basis_index + 1]
        )
        if right_denominator != 0:
            values[:, basis_index] += (
                (knots[basis_index + degree + 1] - x)
                / right_denominator
                * lower_values[:, basis_index + 1]
            )

    return values


def _basis_derivative_matrix_axis(x, knots, degree, derivative_order):
    """Evaluate derivatives of all basis functions for one axis."""
    lower_values = _basis_matrix_axis(
        x,
        knots,
        degree - 1,
        derivative_order=derivative_order - 1,
    )
    n_basis = len(knots) - degree - 1
    values = np.zeros((x.size, n_basis), dtype=float)

    for basis_index in range(n_basis):
        left_denominator = knots[basis_index + degree] - knots[basis_index]
        if left_denominator != 0:
            values[:, basis_index] += (
                degree
                / left_denominator
                * lower_values[:, basis_index]
            )

        right_denominator = (
            knots[basis_index + degree + 1] - knots[basis_index + 1]
        )
        if right_denominator != 0:
            values[:, basis_index] -= (
                degree
                / right_denominator
                * lower_values[:, basis_index + 1]
            )

    return values


def _design_matrix(x, knots, degrees, nu=None, sparse=False):
    """Build the tensor-product design matrix for all dimensions."""
    x = np.asarray(x, dtype=float)
    n_samples, ndim = x.shape
    if nu is None:
        nu = (0,) * ndim

    axis_matrices = [
        _design_matrix_axis(
            x[:, axis],
            knots[axis],
            degrees[axis],
            derivative_order=nu[axis],
        )
        for axis in range(ndim)
    ]

    if sparse:
        return _sparse_design_matrix(axis_matrices, n_samples)

    return _dense_design_matrix(axis_matrices, n_samples)


def _dense_design_matrix(axis_matrices, n_samples):
    """Build a dense tensor-product design matrix."""
    columns = []
    for basis_indices in _basis_index_combinations(axis_matrices):
        column = np.ones(n_samples, dtype=float)

        for axis, basis_index in enumerate(basis_indices):
            column *= axis_matrices[axis][:, basis_index]

        columns.append(column)

    return np.column_stack(columns)


def _sparse_design_matrix(axis_matrices, n_samples):
    """Build a sparse tensor-product design matrix row by row."""
    sparse_module = _require_scipy_sparse()
    basis_counts = [matrix.shape[1] for matrix in axis_matrices]
    column_strides = _column_strides(basis_counts)
    rows = []
    cols = []
    data = []

    for row_index in range(n_samples):
        axis_terms = _nonzero_axis_terms(axis_matrices, row_index)
        if any(len(terms) == 0 for terms in axis_terms):
            continue

        for active_terms in itertools.product(*axis_terms):
            basis_indices = tuple(term[0] for term in active_terms)
            value = np.prod([term[1] for term in active_terms])
            column_index = int(np.dot(basis_indices, column_strides))

            rows.append(row_index)
            cols.append(column_index)
            data.append(value)

    shape = (n_samples, int(np.prod(basis_counts)))
    return sparse_module.csr_matrix((data, (rows, cols)), shape=shape)


def _nonzero_axis_terms(axis_matrices, row_index):
    """Return nonzero basis index/value pairs for one sample row."""
    axis_terms = []
    for matrix in axis_matrices:
        row = matrix[row_index]
        indices = np.nonzero(row)[0]
        axis_terms.append(tuple(zip(indices, row[indices])))

    return axis_terms


def _column_strides(basis_counts):
    """Return C-order flattening strides for tensor-product basis indices."""
    strides = np.empty(len(basis_counts), dtype=int)
    stride = 1
    for axis in range(len(basis_counts) - 1, -1, -1):
        strides[axis] = stride
        stride *= basis_counts[axis]

    return strides


def _basis_index_combinations(axis_matrices):
    """Iterate over tensor-product basis index combinations."""
    basis_counts = [matrix.shape[1] for matrix in axis_matrices]
    return itertools.product(*[range(count) for count in basis_counts])


def _require_scipy_sparse():
    """Import scipy.sparse only when sparse functionality is requested."""
    try:
        from scipy import sparse
    except ImportError as exc:
        raise ImportError("sparse=True requires scipy to be installed.") from exc

    return sparse
