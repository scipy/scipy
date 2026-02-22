"""Test of csgraph public API with int64 index arrays in csr format.

See gh-24629
"""
import numpy as np
import pytest

from scipy.sparse import csr_array
import scipy.sparse.csgraph as csgraph_subpackage
from scipy.sparse.csgraph import (
    connected_components,
    laplacian,
    shortest_path,
    floyd_warshall,
    dijkstra,
    bellman_ford,
    johnson,
    yen,
    breadth_first_order,
    depth_first_order,
    breadth_first_tree,
    depth_first_tree,
    minimum_spanning_tree,
    reverse_cuthill_mckee,
    maximum_flow,
    maximum_bipartite_matching,
    min_weight_full_bipartite_matching,
    structural_rank,
    construct_dist_matrix,
    reconstruct_path,
    csgraph_to_dense,
    csgraph_to_masked,
)


@pytest.fixture()
def A():
    A = csr_array([[0, 1, 2], [2, 1, 0], [0, 1, 0]])
    A.indptr = A.indptr.astype(np.int64)
    A.indices = A.indices.astype(np.int64)
    return A


check_these_functions = [
    # (csgraph_function, args),
    (connected_components, ()),
    (laplacian, ()),
    (floyd_warshall, ()),
    (dijkstra, ()),
    (bellman_ford, ()),
    (johnson, ()),
    (yen, (0, 1, 1)),
    (breadth_first_order, (0,)),
    (depth_first_order, (0,)),
    (breadth_first_tree, (0,)),
    (depth_first_tree, (0,)),
    (minimum_spanning_tree, ()),
    (reverse_cuthill_mckee, ()),
    (maximum_flow, (0, 1)),
    (maximum_bipartite_matching, ()),
    (min_weight_full_bipartite_matching, ()),
    (structural_rank, ()),
    (csgraph_to_dense, ()),
    (csgraph_to_masked, ()),
]

# handled individually outside of check_these_function
handled = ["shortest_path", "construct_dist_matrix", "reconstruct_path",
           "csgraph_from_dense", "csgraph_from_masked", "csgraph_masked_from_dense"]

# look for functions that should be tested
functions = [fname for fname in dir(csgraph_subpackage)
             if fname[0] != "_"
             if not fname.endswith("Error")
             if not fname.startswith("test")
             if fname not in handled]

checked_fs = [x[0].__name__ for x in check_these_functions]
for f in functions:
    if f not in checked_fs:
        raise ValueError(f"function {f} in csgraph is not tested for index arrays")


@pytest.mark.parametrize("csgraph_function, args", check_these_functions)
def test_smoke_from_int64_index_arrays(csgraph_function, args, A):
    csgraph_function(A, *args)


def test_smoke_from_int64_construct_dist_matrix(A):
    _, preds = shortest_path(A, return_predecessors=True)
    construct_dist_matrix(A, preds)


def test_smoke_from_int64_shortest_path(A):
    shortest_path(A, method='auto')


def test_smoke_from_int64_reconstruct_path(A):
    pred = np.array([-9999, 0, 1])
    reconstruct_path(A, pred)
