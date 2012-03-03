import numpy as np
from numpy.testing import assert_array_almost_equal
from scipy.sparse.csgraph import \
    shortest_path, dijkstra, floyd_warshall, bellman_ford, construct_dist_matrix


def floyd_warshall_slow(graph, directed=False):
    N = graph.shape[0]

    #set nonzero entries to infinity
    graph[np.where(graph == 0)] = np.inf

    #set diagonal to zero
    graph.flat[::N + 1] = 0

    if not directed:
        graph = np.minimum(graph, graph.T)

    for k in range(N):
        for i in range(N):
            for j in range(N):
                graph[i, j] = min(graph[i, j], graph[i, k] + graph[k, j])

    return graph


def generate_graph(N=20):
    #sparse grid of distances
    dist_matrix = np.random.random((N, N))

    #make graph sparse
    i = (np.random.randint(N, size=N * N / 2),
         np.random.randint(N, size=N * N / 2))
    dist_matrix[i] = 0

    #set diagonal to zero
    dist_matrix.flat[::N + 1] = 0

    return dist_matrix


def test_floyd_warshall():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_FW = shortest_path(dist_matrix, 'FW', directed)
        graph_py = floyd_warshall_slow(dist_matrix.copy(), directed)

        assert_array_almost_equal(graph_FW, graph_py)


def test_dijkstra():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_D = shortest_path(dist_matrix, 'D', directed)
        graph_FW = shortest_path(dist_matrix, 'FW', directed)

        assert_array_almost_equal(graph_D, graph_FW)


def test_dijkstra_indices():
    dist_matrix = generate_graph(20)
    indices = np.arange(20, dtype=int)
    np.random.seed(0)
    np.random.shuffle(indices)
    indices = indices[:6]

    for directed in (True, False):
        graph_FW = shortest_path(dist_matrix, 'FW', directed)
        for indshape in [(6,), (6, 1), (2, 3)]:
            outshape = indshape + (20,)
            graph_D = dijkstra(dist_matrix, directed,
                               indices=indices.reshape(indshape))
            assert_array_almost_equal(graph_D,
                                      graph_FW[indices].reshape(outshape))


def test_johnson():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_J = shortest_path(dist_matrix, 'J', directed)
        graph_FW = shortest_path(dist_matrix, 'FW', directed)

        assert_array_almost_equal(graph_J, graph_FW)


def test_johnson_indices():
    dist_matrix = generate_graph(20)
    indices = np.arange(20, dtype=int)
    np.random.seed(0)
    np.random.shuffle(indices)
    indices = indices[:6]

    for directed in (True, False):
        graph_FW = shortest_path(dist_matrix, 'FW', directed)
        for indshape in [(6,), (6, 1), (2, 3)]:
            outshape = indshape + (20,)
            graph_J = dijkstra(dist_matrix, directed,
                               indices=indices.reshape(indshape))
            assert_array_almost_equal(graph_J,
                                      graph_FW[indices].reshape(outshape))


def test_bellman_ford():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_BF = shortest_path(dist_matrix, 'BF', directed)
        graph_py = floyd_warshall_slow(dist_matrix.copy(), directed)

        assert_array_almost_equal(graph_BF, graph_py)


def test_bellman_ford_indices():
    dist_matrix = generate_graph(20)
    indices = np.arange(20, dtype=int)
    np.random.seed(0)
    np.random.shuffle(indices)
    indices = indices[:6]

    for directed in (True, False):
        graph_FW = shortest_path(dist_matrix, 'FW', directed)
        for indshape in [(6,), (6, 1), (2, 3)]:
            outshape = indshape + (20,)
            graph_BF = bellman_ford(dist_matrix, directed,
                                    indices=indices.reshape(indshape))
            assert_array_almost_equal(graph_BF,
                                      graph_FW[indices].reshape(outshape))


def test_predecessors():
    csgraph = generate_graph(20)

    for directed in (True, False):
        dist_D, pred_D = shortest_path(csgraph, 'D', directed,
                                       return_predecessors=True)
        dist_FW, pred_FW = shortest_path(csgraph, 'FW', directed,
                                         return_predecessors=True)

        assert_array_almost_equal(dist_D, dist_FW)
        assert_array_almost_equal(pred_D, pred_FW)


def test_construct_shortest_path():
    csgraph = generate_graph(5)

    for directed in (True, False):
        dist, pred = shortest_path(csgraph,
                                   directed=directed,
                                   overwrite=False,
                                   return_predecessors=True)
        dist2 = construct_dist_matrix(csgraph, pred, directed=directed)

        assert_array_almost_equal(dist, dist2)

def test_unweighted_path():
    csgraph = generate_graph(20)
    csgraph_ones = np.ones(csgraph.shape)
    csgraph_ones[csgraph == 0] = 0

    for directed in (True, False):
        for method in ('FW', 'D'):
            D1 = shortest_path(csgraph, method=method, directed=directed,
                               unweighted=True)
            D2 = shortest_path(csgraph_ones, method=method, directed=directed)
            assert_array_almost_equal(D1, D2)
        


if __name__ == '__main__':
    import nose
    nose.runmodule()
