import numpy as np
from .common import Benchmark

try:
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial.distance import cdist
except ImportError:
    pass


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
    param_names = ['shapes', 'cost_type']
    params = [
        [(i, i) for i in sizes] + [(i, 2 * i) for i in sizes],
        ['uniform', 'spatial', 'logarithmic', 'integer', 'binary']
    ]

    def setup(self, shape, cost_type):

        cost_func = {'uniform': random_uniform,
                     'spatial': random_spatial,
                     'logarithmic': random_logarithmic,
                     'integer': random_integer,
                     'binary': random_binary}[cost_type]

        self.cost_matrix = cost_func(shape)

    def time_evaluation(self, *args):
        linear_sum_assignment(self.cost_matrix)
