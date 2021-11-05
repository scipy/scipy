import numpy as np
from numpy.testing import assert_equal
from pytest import raises as assert_raises
from scipy.sparse.csgraph import (
    csgraph_from_adjacency_list,
    csgraph_to_adjacency_list,
    csgraph_from_dense,
)


def test_csgraph_weighted_adjacency_list():
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

    unweighted_graph_adjacency_list = {
        0: {1: 1, 2: 1},
        1: {3: 1},
        2: {3: 1},
        3: {},
    }

    unweighted_graph_matrix = np.array([
        [0, 1, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
    ])

    # test csgraph_from_adjacency_list
    for weighted in [True, False]:
        expected = graph_matrix if weighted else unweighted_graph_matrix

        assert_equal(
            csgraph_from_adjacency_list(
                graph_adjacency_list, weighted=weighted
            ).toarray(),
            expected,
        )

        assert_equal(
            csgraph_from_adjacency_list(
                graph_adjacency_list, weighted=weighted
            ).toarray(),
            csgraph_from_dense(expected).toarray(),
        )

    # test csgraph_to_adjacency_list
    for weighted in [True, False]:
        expected = graph_adjacency_list if weighted else unweighted_graph_adjacency_list

        assert_equal(
            csgraph_to_adjacency_list(
                csgraph_from_adjacency_list(
                    graph_adjacency_list, weighted=weighted
                )
            ),
            expected,
        )

    # test csgraph_to_adjacency_list with csgraph_from_dense
    assert_equal(
        csgraph_to_adjacency_list(
            csgraph_from_dense(graph_matrix)
        ),
        graph_adjacency_list,
    )

    assert_equal(
        csgraph_to_adjacency_list(
            csgraph_from_dense(unweighted_graph_matrix)
        ),
        unweighted_graph_adjacency_list,
    )

    assert_raises(
        ValueError,
        csgraph_from_adjacency_list,
        {
            0: [1, 2],
            1: [3],
            2: [3],
            3: [],
        },
        weighted=True,
    )
