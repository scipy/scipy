import numpy as np

from scipy.linalg import norm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import Delaunay
from scipy.spatial.distance import squareform, pdist


def _delaunay_to_euclidean_graph(tri):
    """TODO"""

    npoints, ndim = tri.points.shape

    # edge to weight dict
    edge_to_weight = {}
    for simplex in tri.simplices:
        for i in range(-1, ndim):
            node_from, node_to = sorted(
                [simplex[i], simplex[i+1]]
            )  # for undirected graph
            if (node_from, node_to) in edge_to_weight:
                continue

            weight = norm(tri.points[node_from] - tri.points[node_to])
            edge_to_weight[(node_from, node_to)] = weight

    # edge_to_weight to csgraph (csr_matrix)
    n_edges = len(edge_to_weight)

    row = np.zeros(n_edges, dtype=int)
    col = np.zeros(n_edges, dtype=int)
    data = np.zeros(n_edges, dtype=float)

    for i, ((node_from, node_to), weight) in enumerate(edge_to_weight.items()):
        row[i] = node_from
        col[i] = node_to
        data[i] = weight

    return csr_matrix((data, (row, col)), shape=(npoints, npoints))


def euclidean_minimum_spanning_tree(points, method="delaunay"):
    """TODO"""

    if method == "pairwise":
        graph = squareform(pdist(points))

        # undirected graph
        inds = np.tril_indices(points.shape[0])
        graph[inds] = 0

        graph = csr_matrix(graph)

    elif method == "delaunay":
        tri = Delaunay(points)
        graph = _delaunay_to_euclidean_graph(tri)

    else:
        raise ValueError("unrecognized method '%s'" % method)

    return minimum_spanning_tree(graph)
