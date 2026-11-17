"""Regression tests for scipy/scipy#26549 and matrix-free narrow products."""

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from scipy import sparse
from scipy.sparse.csgraph import laplacian


def apply(operator, v):
    return operator(v) if callable(operator) else operator @ v


@pytest.mark.parametrize("factory", [np.array, sparse.csr_array, sparse.coo_array])
@pytest.mark.parametrize("form", ["array", "function", "lo"])
@pytest.mark.parametrize("symmetrized", [False, True])
@pytest.mark.parametrize("copy", [False, True])
def test_requested_wide_dtype_preserves_degrees(factory, form, symmetrized, copy):
    graph = factory(np.array([[0, 100, 100], [100, 0, 0], [100, 0, 0]], dtype=np.int8))
    original = graph.copy()
    expected = np.array([[200, -100, -100], [-100, 100, 0], [-100, 0, 100]])
    if symmetrized:
        expected *= 2
    operator, degree = laplacian(
        graph, dtype=np.int64, form=form,
        symmetrized=symmetrized, copy=copy, return_diag=True,
    )
    assert_array_equal(apply(operator, np.eye(3, dtype=np.int8)), expected)
    assert_array_equal(apply(operator, np.ones((3, 1), dtype=np.int8)), np.zeros((3, 1)))
    assert_array_equal(degree, np.diag(expected))
    assert_array_equal(
        graph.toarray() if sparse.issparse(graph) else graph,
        original.toarray() if sparse.issparse(original) else original,
    )


@pytest.mark.parametrize("factory", [np.array, sparse.csr_array, sparse.coo_array])
@pytest.mark.parametrize("form", ["array", "function", "lo"])
@pytest.mark.parametrize("copy", [False, True])
def test_wide_dtype_product_without_degree_overflow(factory, form, copy):
    graph = factory(np.array([[0, 100], [100, 0]], dtype=np.int8))
    operator = laplacian(graph, dtype=np.int64, form=form, copy=copy)
    vector = np.array([[0], [65]], dtype=np.int8)
    assert_array_equal(apply(operator, vector), [[-6500], [6500]])


@pytest.mark.parametrize("dtype", [np.uint8, np.uint16, np.uint32, np.uint64])
@pytest.mark.parametrize("form", ["array", "function", "lo"])
def test_normalized_unsigned_graph(dtype, form):
    graph = np.array([[0, 1], [1, 0]], dtype=dtype)
    operator = laplacian(graph, normed=True, form=form)
    assert_allclose(apply(operator, np.eye(2)), [[1, -1], [-1, 1]])


def test_copy_false_preserves_in_place_path_when_no_promotion():
    graph = np.array([[0, 4], [4, 0]], dtype=np.int64)
    result = laplacian(graph, dtype=np.int64, copy=False)
    assert np.shares_memory(result, graph)
    assert_array_equal(graph, [[4, -4], [-4, 4]])
