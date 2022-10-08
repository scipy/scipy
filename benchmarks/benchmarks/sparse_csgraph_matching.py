import numpy as np
import scipy.sparse
from scipy.spatial.distance import cdist

from .common import Benchmark, safe_import


with safe_import():
    from scipy.sparse.csgraph import maximum_bipartite_matching,\
        min_weight_full_bipartite_matching


class MaximumBipartiteMatching(Benchmark):
    params = [[5000, 7500, 10000], [0.0001, 0.0005, 0.001]]
    param_names = ['n', 'density']

    def setup(self, n, density):
        # Create random sparse matrices. Note that we could use
        # scipy.sparse.rand for this purpose, but simply using np.random and
        # disregarding duplicates is quite a bit faster.
        rng = np.random.default_rng(42)
        d = rng.integers(0, n, size=(int(n*n*density), 2))
        graph = scipy.sparse.csr_matrix((np.ones(len(d)), (d[:, 0], d[:, 1])),
                                        shape=(n, n))
        self.graph = graph

    def time_maximum_bipartite_matching(self, n, density):
        maximum_bipartite_matching(self.graph)


# For benchmarking min_weight_full_bipartite_matching, we rely on some of
# the classes defined in Burkard, Dell'Amico, Martello -- Assignment Problems,
# 2009, Section 4.10.1.
def random_uniform(shape, rng):
    return scipy.sparse.csr_matrix(rng.uniform(1, 100, shape))


def random_uniform_sparse(shape, rng):
    return scipy.sparse.random(shape[0], shape[1], density=0.1, format='csr', random_state=rng)


def random_uniform_integer(shape, rng):
    return scipy.sparse.csr_matrix(rng.integers(1, 1000, shape))


def random_geometric(shape, rng):
    P = rng.integers(1, 1000, size=(shape[0], 2))
    Q = rng.integers(1, 1000, size=(shape[1], 2))
    return scipy.sparse.csr_matrix(cdist(P, Q, 'sqeuclidean'))


def random_two_cost(shape, rng):
    return scipy.sparse.csr_matrix(rng.choice((1, 1000000), shape))


def machol_wien(shape, rng):
    # Machol--Wien instances being harder than the other examples, we cut
    # down the size of the instance by 5.
    return scipy.sparse.csr_matrix(
        np.outer(np.arange(shape[0]//5) + 1, np.arange(shape[1]//5) + 1))


class MinWeightFullBipartiteMatching(Benchmark):

    sizes = range(100, 401, 100)
    param_names = ['shapes', 'input_type']
    params = [
        [(i, i) for i in sizes] + [(i, 2 * i) for i in sizes],
        ['random_uniform', 'random_uniform_sparse', 'random_uniform_integer',
         'random_geometric', 'random_two_cost', 'machol_wien']
    ]

    def setup(self, shape, input_type):
        rng = np.random.default_rng(42)
        input_func = {'random_uniform': random_uniform,
                      'random_uniform_sparse': random_uniform_sparse,
                      'random_uniform_integer': random_uniform_integer,
                      'random_geometric': random_geometric,
                      'random_two_cost': random_two_cost,
                      'machol_wien': machol_wien}[input_type]

        self.biadjacency_matrix = input_func(shape, rng)

    def time_evaluation(self, *args):
        min_weight_full_bipartite_matching(self.biadjacency_matrix)
