import numpy as np
from numpy.testing import assert_equal
from scipy.sparse.csgraph import (
    csgraph_from_adjacency_list,
    csgraph_to_adjacency_list,
    csgraph_from_dense,
)


def test_csgraph_adjacency_list():
    graph_adjacency_list = {
        0: {1: 1, 2: 2},
        1: {3: 1},
        2: {3: 3},
        3: {},
    }

    graph_matrix = np.array([
        [0, 1, 2, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 3],
        [0, 0, 0, 0],
    ])

    assert_equal(
        csgraph_from_adjacency_list(graph_adjacency_list).toarray(),
        graph_matrix,
    )
    assert_equal(
        csgraph_from_adjacency_list(graph_adjacency_list).toarray(),
        csgraph_from_dense(graph_matrix).toarray(),
    )
    assert_equal(
        csgraph_to_adjacency_list(
            csgraph_from_adjacency_list(graph_adjacency_list)
        ),
        graph_adjacency_list,
    )
    assert_equal(
        csgraph_to_adjacency_list(
            csgraph_from_dense(graph_matrix)
        ),
        graph_adjacency_list,
    )
