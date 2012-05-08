import numpy as np
from numpy.testing import assert_array_almost_equal, assert_raises, TestCase
from scipy.sparse.csgraph import \
    shortest_path, dijkstra, floyd_warshall, johnson,\
    bellman_ford, construct_dist_matrix, NegativeCycleError


directed_G = np.array([[0, 3, 3, 0, 0],
                       [0, 0, 0, 2, 4],
                       [0, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0],
                       [2, 0, 0, 2, 0]], dtype=float)

undirected_G = np.array([[0, 3, 3, 1, 2],
                         [3, 0, 0, 2, 4],
                         [3, 0, 0, 0, 0],
                         [1, 2, 0, 0, 2],
                         [2, 4, 0, 2, 0]], dtype=float)

unweighted_G = (directed_G > 0).astype(float)

directed_SP = [[0, 3, 3, 5, 7],
               [3, 0, 6, 2, 4],
               [np.inf, np.inf, 0, np.inf, np.inf],
               [1, 4, 4, 0, 8],
               [2, 5, 5, 2, 0]]

directed_pred = np.array([[-9999,     0,     0,     1,     1],
                          [    3, -9999,     0,     1,     1],
                          [-9999, -9999, -9999, -9999, -9999],
                          [    3,     0,     0, -9999,     1],
                          [    4,     0,     0,     4, -9999]], dtype=float)

undirected_SP = np.array([[0, 3, 3, 1, 2],
                          [3, 0, 6, 2, 4],
                          [3, 6, 0, 4, 5],
                          [1, 2, 4, 0, 2],
                          [2, 4, 5, 2, 0]], dtype=float)

undirected_pred = np.array([[-9999,     0,     0,     0,     0],
                            [    1, -9999,     0,     1,     1],
                            [    2,     0, -9999,     0,     0],
                            [    3,     3,     0, -9999,     3],
                            [    4,     4,     0,     4, -9999]], dtype=float)


methods = ['auto', 'FW', 'D', 'BF', 'J']

def test_directed():
    for method in methods:
        SP = shortest_path(directed_G, method=method, directed=True,
                           overwrite=False)
        yield (assert_array_almost_equal, SP, directed_SP)

def test_undirected():
    for method in methods:
        SP1 = shortest_path(directed_G, method=method, directed=False,
                            overwrite=False)
        SP2 = shortest_path(undirected_G, method=method, directed=True,
                            overwrite=False)

        yield (assert_array_almost_equal, SP1, undirected_SP)
        yield (assert_array_almost_equal, SP2, undirected_SP)

def test_shortest_path_indices():
    indices = np.arange(4)
    for indshape in [(4,), (4, 1), (2, 2)]:
        for func in (dijkstra, johnson, bellman_ford):
            outshape = indshape + (5,)
            SP = func(directed_G, directed=False,
                      indices=indices.reshape(indshape))
            yield (assert_array_almost_equal, SP,
                   undirected_SP[indices].reshape(outshape))

def test_predecessors():
    SP_res = {True: directed_SP,
              False: undirected_SP}
    pred_res = {True: directed_pred,
                False: undirected_pred}
    
    for method in methods:
        for directed in True, False:
            SP, pred = shortest_path(directed_G, method, directed=directed,
                                     overwrite=False,
                                     return_predecessors=True)

            yield (assert_array_almost_equal,
                   SP, SP_res[directed])
            yield (assert_array_almost_equal,
                   pred, pred_res[directed])

def test_construct_shortest_path():
    SP_res = {True: directed_SP,
              False: undirected_SP}
    for method in methods:
        for directed in (True, False):
            SP1, pred = shortest_path(directed_G,
                                      directed=directed,
                                      overwrite=False,
                                      return_predecessors=True)
            SP2 = construct_dist_matrix(directed_G, pred, directed=directed)

            yield (assert_array_almost_equal, SP1, SP2)

def test_unweighted_path():
    for method in methods:
        for directed in (True, False):
            SP1 = shortest_path(directed_G,
                                directed=directed,
                                overwrite=False,
                                unweighted=True)
            SP2 = shortest_path(unweighted_G,
                                directed=directed,
                                overwrite=False,
                                unweighted=False)
    
            yield (assert_array_almost_equal, SP1, SP2)


def test_negative_cycles():
    # create a small graph with a negative cycle
    graph = np.ones([5, 5])
    graph.flat[::6] = 0
    graph[1, 2] = -2
    for method in ['FW', 'J', 'BF']:
        for directed in True, False:
            yield (assert_raises, NegativeCycleError, shortest_path,
                   graph, method, directed)

def test_masked_input():
    G = np.ma.masked_equal(directed_G, 0)
    for method in methods:
        SP = shortest_path(directed_G, method=method, directed=True,
                           overwrite=False)
        yield (assert_array_almost_equal, SP, directed_SP)
    

if __name__ == '__main__':
    import nose
    nose.runmodule()
