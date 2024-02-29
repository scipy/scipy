from concurrent.futures import ThreadPoolExecutor, wait

import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial.distance import cdist


def random_uniform(shape):
    return np.random.uniform(-20, 20, shape)


def random_logarithmic(shape):
    return 10**np.random.uniform(-20, 20, shape)


def random_integer(shape):
    return np.random.randint(-1000, 1000, shape)


def random_binary(shape):
    return np.random.randint(0, 2, shape)


def random_spatial(shape):
    P = np.random.uniform(-1, 1, size=(shape[0], 2))
    Q = np.random.uniform(-1, 1, size=(shape[1], 2))
    return cdist(P, Q, 'sqeuclidean')


class LinearAssignment(Benchmark):

    sizes = range(100, 401, 100)
    shapes = [(i, i) for i in sizes]
    shapes.extend([(i, 2 * i) for i in sizes])
    shapes.extend([(2 * i, i) for i in sizes])
    cost_types = ['uniform', 'spatial', 'logarithmic', 'integer', 'binary']
    param_names = ['shape', 'cost_type']
    params = [shapes, cost_types]

    def setup(self, shape, cost_type):

        cost_func = {'uniform': random_uniform,
                     'spatial': random_spatial,
                     'logarithmic': random_logarithmic,
                     'integer': random_integer,
                     'binary': random_binary}[cost_type]

        self.cost_matrix = cost_func(shape)

    def time_evaluation(self, *args):
        linear_sum_assignment(self.cost_matrix)


class ParallelLinearAssignment(Benchmark):
    shape = (100, 100)
    param_names = ['threads']
    params = [[1, 2, 4]]

    def setup(self, threads):
        self.cost_matrices = [random_uniform(self.shape) for _ in range(20)]

    def time_evaluation(self, threads):
        with ThreadPoolExecutor(max_workers=threads) as pool:
            wait({pool.submit(linear_sum_assignment, cost_matrix)
                  for cost_matrix in self.cost_matrices})
