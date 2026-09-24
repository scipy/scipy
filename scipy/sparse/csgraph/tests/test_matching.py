from itertools import product

import numpy as np
from numpy.testing import assert_array_equal, assert_equal
import pytest

from scipy.sparse import csr_array, csc_array, coo_array, diags_array
from scipy.sparse.csgraph import (
    maximum_bipartite_matching, min_weight_full_bipartite_matching,
    structural_rank
)


def test_maximum_bipartite_matching_raises_on_dense_input():
    with pytest.raises(TypeError):
        graph = np.array([[0, 1], [0, 0]])
        maximum_bipartite_matching(graph)


def test_maximum_bipartite_matching_empty_graph():
    graph = csr_array((0, 0))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    expected_matching = np.array([])
    assert_array_equal(expected_matching, x)
    assert_array_equal(expected_matching, y)


def test_maximum_bipartite_matching_empty_left_partition():
    graph = csr_array((2, 0))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    assert_array_equal(np.array([]), x)
    assert_array_equal(np.array([-1, -1]), y)


def test_maximum_bipartite_matching_empty_right_partition():
    graph = csr_array((0, 3))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    assert_array_equal(np.array([-1, -1, -1]), x)
    assert_array_equal(np.array([]), y)


def test_maximum_bipartite_matching_graph_with_no_edges():
    graph = csr_array((2, 2))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    assert_array_equal(np.array([-1, -1]), x)
    assert_array_equal(np.array([-1, -1]), y)


def test_maximum_bipartite_matching_graph_that_causes_augmentation():
    # In this graph, column 1 is initially assigned to row 1, but it should be
    # reassigned to make room for row 2.
    graph = csr_array([[1, 1], [1, 0]])
    x = maximum_bipartite_matching(graph, perm_type='column')
    y = maximum_bipartite_matching(graph, perm_type='row')
    expected_matching = np.array([1, 0])
    assert_array_equal(expected_matching, x)
    assert_array_equal(expected_matching, y)


def test_maximum_bipartite_matching_graph_with_more_rows_than_columns():
    graph = csr_array([[1, 1], [1, 0], [0, 1]])
    x = maximum_bipartite_matching(graph, perm_type='column')
    y = maximum_bipartite_matching(graph, perm_type='row')
    assert_array_equal(np.array([0, -1, 1]), x)
    assert_array_equal(np.array([0, 2]), y)


def test_maximum_bipartite_matching_graph_with_more_columns_than_rows():
    graph = csr_array([[1, 1, 0], [0, 0, 1]])
    x = maximum_bipartite_matching(graph, perm_type='column')
    y = maximum_bipartite_matching(graph, perm_type='row')
    assert_array_equal(np.array([0, 2]), x)
    assert_array_equal(np.array([0, -1, 1]), y)


def test_maximum_bipartite_matching_explicit_zeros_count_as_edges():
    data = [0, 0]
    indices = [1, 0]
    indptr = [0, 1, 2]
    graph = csr_array((data, indices, indptr), shape=(2, 2))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    expected_matching = np.array([1, 0])
    assert_array_equal(expected_matching, x)
    assert_array_equal(expected_matching, y)


def test_maximum_bipartite_matching_feasibility_of_result():
    # This is a regression test for GitHub issue #11458
    data = np.ones(50, dtype=int)
    indices = [11, 12, 19, 22, 23, 5, 22, 3, 8, 10, 5, 6, 11, 12, 13, 5, 13,
               14, 20, 22, 3, 15, 3, 13, 14, 11, 12, 19, 22, 23, 5, 22, 3, 8,
               10, 5, 6, 11, 12, 13, 5, 13, 14, 20, 22, 3, 15, 3, 13, 14]
    indptr = [0, 5, 7, 10, 10, 15, 20, 22, 22, 23, 25, 30, 32, 35, 35, 40, 45,
              47, 47, 48, 50]
    graph = csr_array((data, indices, indptr), shape=(20, 25))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    assert (x != -1).sum() == 13
    assert (y != -1).sum() == 13
    # Ensure that each element of the matching is in fact an edge in the graph.
    for u, v in zip(range(graph.shape[0]), y):
        if v != -1:
            assert graph[u, v]
    for u, v in zip(x, range(graph.shape[1])):
        if u != -1:
            assert graph[u, v]


def test_maximum_bipartite_matching_duplicate_entries_in_csr_input():
    # Regression test for GitHub issue #26160: a CSR input that stores the
    # same entry more than once (so has_canonical_format is False) used to
    # push the same row onto the fixed-size DFS stack repeatedly, overflowing
    # it and corrupting memory. Rows are {0, 1} and {0}, with column 0 stored
    # K extra times in row 1.
    K = 100_000
    data = np.ones(2 + K)
    indices = [0, 1] + [0] * K
    indptr = [0, 2, 2 + K]
    graph = csr_array((data, indices, indptr), shape=(2, 2))
    assert not graph.has_canonical_format
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    expected_matching = np.array([1, 0])
    assert_array_equal(x, expected_matching)
    assert_array_equal(y, expected_matching)


def test_structural_rank_duplicate_entries_in_csr_input():
    # Structural rank delegates to maximum_bipartite_matching, so it must be
    # guarded against duplicate CSR entries the same way (gh-26160).
    K = 100_000
    data = np.ones(2 + K)
    indices = [0, 1] + [0] * K
    indptr = [0, 2, 2 + K]
    graph = csr_array((data, indices, indptr), shape=(2, 2))
    assert_equal(structural_rank(graph), 2)


def test_min_weight_full_bipartite_matching_duplicate_entries_in_csr_input():
    # Regression test for GitHub issue #26160: min_weight_full_bipartite_matching
    # used to overflow fixed-size stacks in _hopcroft_karp and _lapjvsp when the
    # (weighted) input stored the same entry more than once.
    K = 100_000
    data = np.ones(2 + K)
    indices = [0, 1] + [0] * K
    indptr = [0, 2, 2 + K]
    graph = csr_array((data, indices, indptr), shape=(2, 2))
    src, dst = min_weight_full_bipartite_matching(graph)
    # The graph is 2 x 2 with an edge from each row, so a full matching must
    # exist and cover both rows and both columns.
    assert_array_equal(src, np.array([0, 1]))
    assert_array_equal(dst, np.array([1, 0]))
    for u, v in zip(src, dst):
        assert graph[u, v] != 0


def test_duplicate_entries_in_csc_input():
    # gh-26160: the duplicate-entry overflow is also reachable through a CSC
    # input, because converting CSC to CSR preserves duplicates and the
    # column-oriented twin of the CSR test graph has column 0 connected to
    # both rows (row 1 K times) and column 1 to row 0.
    K = 100_000
    data = np.ones(2 + K)
    indices = [0] + [1] * K + [0]
    indptr = [0, 1 + K, 2 + K]
    graph = csc_array((data, indices, indptr), shape=(2, 2))
    assert not graph.has_canonical_format
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    expected_matching = np.array([1, 0])
    assert_array_equal(x, expected_matching)
    assert_array_equal(y, expected_matching)
    assert_equal(structural_rank(graph), 2)
    src, dst = min_weight_full_bipartite_matching(graph)
    assert_array_equal(src, np.array([0, 1]))
    assert_array_equal(dst, np.array([1, 0]))


def test_matching_input_duplicate_entries_not_mutated():
    # gh-26160: duplicate entries must be merged on an internal copy so that
    # the caller's matrix is left untouched. tocsr() returns a CSR input
    # unchanged and sum_duplicates() mutates in place.
    K = 100_000
    data = np.ones(2 + K)
    indices = [0, 1] + [0] * K
    indptr = [0, 2, 2 + K]
    for func in (maximum_bipartite_matching, min_weight_full_bipartite_matching):
        graph = csr_array((data, indices, indptr), shape=(2, 2))
        data_before = graph.data.copy()
        indices_before = graph.indices.copy()
        indptr_before = graph.indptr.copy()
        canonical_before = graph.has_canonical_format
        if func is maximum_bipartite_matching:
            func(graph, perm_type='row')
        else:
            func(graph)
        assert graph.has_canonical_format == canonical_before
        assert_array_equal(graph.data, data_before)
        assert_array_equal(graph.indices, indices_before)
        assert_array_equal(graph.indptr, indptr_before)


def test_matching_large_random_graph_with_one_edge_incident_to_each_vertex():
    np.random.seed(42)
    A = diags_array(np.ones(25), offsets=0, format='csr')
    rand_perm = np.random.permutation(25)
    rand_perm2 = np.random.permutation(25)

    Rrow = np.arange(25)
    Rcol = rand_perm
    Rdata = np.ones(25, dtype=int)
    Rmat = csr_array((Rdata, (Rrow, Rcol)))

    Crow = rand_perm2
    Ccol = np.arange(25)
    Cdata = np.ones(25, dtype=int)
    Cmat = csr_array((Cdata, (Crow, Ccol)))
    # Randomly permute identity matrix
    B = Rmat @ A @ Cmat

    # Row permute
    perm = maximum_bipartite_matching(B, perm_type='row')
    Rrow = np.arange(25)
    Rcol = perm
    Rdata = np.ones(25, dtype=int)
    Rmat = csr_array((Rdata, (Rrow, Rcol)))
    C1 = Rmat @ B

    # Column permute
    perm2 = maximum_bipartite_matching(B, perm_type='column')
    Crow = perm2
    Ccol = np.arange(25)
    Cdata = np.ones(25, dtype=int)
    Cmat = csr_array((Cdata, (Crow, Ccol)))
    C2 = B @ Cmat

    # Should get identity matrix back
    assert_equal(any(C1.diagonal() == 0), False)
    assert_equal(any(C2.diagonal() == 0), False)


@pytest.mark.parametrize('num_rows,num_cols', [(0, 0), (2, 0), (0, 3)])
def test_min_weight_full_matching_trivial_graph(num_rows, num_cols):
    biadjacency = csr_array((num_cols, num_rows))
    biadjacency1 = coo_array((num_cols, num_rows))

    row_ind, col_ind = min_weight_full_bipartite_matching(biadjacency)
    assert len(row_ind) == 0
    assert len(col_ind) == 0

    row_ind1, col_ind1 = min_weight_full_bipartite_matching(biadjacency1)
    assert len(row_ind1) == 0
    assert len(col_ind1) == 0


@pytest.mark.parametrize('biadjacency',
                         [
                            [[1, 1, 1], [1, 0, 0], [1, 0, 0]],
                            [[1, 1, 1], [0, 0, 1], [0, 0, 1]],
                            [[1, 0, 0, 1], [1, 1, 0, 1], [0, 0, 0, 0]],
                            [[1, 0, 0], [2, 0, 0]],
                            [[0, 1, 0], [0, 2, 0]],
                            [[1, 0], [2, 0], [5, 0]]
                         ])
def test_min_weight_full_matching_infeasible_problems(biadjacency):
    with pytest.raises(ValueError):
        min_weight_full_bipartite_matching(csr_array(biadjacency))
    with pytest.raises(ValueError):
        min_weight_full_bipartite_matching(coo_array(biadjacency))


def test_min_weight_full_matching_large_infeasible():
    # Regression test for GitHub issue #17269
    a = np.asarray([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001],
        [0.0, 0.11687445, 0.0, 0.0, 0.01319788, 0.07509257, 0.0,
         0.0, 0.0, 0.74228317, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.81087935, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.8408466, 0.0, 0.0, 0.0, 0.0, 0.01194389,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.82994211, 0.0, 0.0, 0.0, 0.11468516, 0.0, 0.0, 0.0,
         0.11173505, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0],
        [0.18796507, 0.0, 0.04002318, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75883335,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.71545464, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02748488,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.78470564, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.14829198,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.10870609, 0.0, 0.0, 0.0, 0.8918677, 0.0, 0.0, 0.0, 0.06306644,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.63844085, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7442354, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.09850549, 0.0, 0.0, 0.18638258,
         0.2769244, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.73182464, 0.0, 0.0, 0.46443561,
         0.38589284, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.29510278, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.09666032, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ])
    with pytest.raises(ValueError, match='no full matching exists'):
        min_weight_full_bipartite_matching(csr_array(a))
    with pytest.raises(ValueError, match='no full matching exists'):
        min_weight_full_bipartite_matching(coo_array(a))


def test_explicit_zero_causes_warning():
    biadjacency = csr_array(((2, 0, 3), (0, 1, 1), (0, 2, 3)))
    with pytest.warns(UserWarning):
        min_weight_full_bipartite_matching(biadjacency)
    with pytest.warns(UserWarning):
        min_weight_full_bipartite_matching(biadjacency.tocoo())


# General test for linear sum assignment solvers to make it possible to rely
# on the same tests for scipy.optimize.linear_sum_assignment.
def linear_sum_assignment_assertions(
    solver, array_type, sign, test_case
):
    cost_matrix, expected_cost = test_case
    maximize = sign == -1
    cost_matrix = sign * array_type(cost_matrix)
    expected_cost = sign * np.array(expected_cost)

    row_ind, col_ind = solver(cost_matrix, maximize=maximize)
    assert_array_equal(row_ind, np.sort(row_ind))
    assert_array_equal(expected_cost,
                       np.array(cost_matrix[row_ind, col_ind]).flatten())

    cost_matrix = cost_matrix.T
    row_ind, col_ind = solver(cost_matrix, maximize=maximize)
    assert_array_equal(row_ind, np.sort(row_ind))
    assert_array_equal(np.sort(expected_cost),
                       np.sort(np.array(
                           cost_matrix[row_ind, col_ind])).flatten())


linear_sum_assignment_test_cases = list(product(
    [-1, 1],
    [
        # Square
        ([[400, 150, 400],
          [400, 450, 600],
          [300, 225, 300]],
         [150, 400, 300]),

        # Rectangular variant
        ([[400, 150, 400, 1],
          [400, 450, 600, 2],
          [300, 225, 300, 3]],
         [150, 2, 300]),

        ([[10, 10, 8],
          [9, 8, 1],
          [9, 7, 4]],
         [10, 1, 7]),

        # Square
        ([[10, 10, 8, 11],
          [9, 8, 1, 1],
          [9, 7, 4, 10]],
         [10, 1, 4]),

        # Rectangular variant
        ([[10, float("inf"), float("inf")],
          [float("inf"), float("inf"), 1],
          [float("inf"), 7, float("inf")]],
         [10, 1, 7])
    ]))


@pytest.mark.parametrize('sign,test_case', linear_sum_assignment_test_cases)
def test_min_weight_full_matching_small_inputs(sign, test_case):
    linear_sum_assignment_assertions(
        min_weight_full_bipartite_matching, csr_array, sign, test_case)
