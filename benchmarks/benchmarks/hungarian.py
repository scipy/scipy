
try:
    from numpy.random import rand
    import numpy.random as nprand
    from scipy.optimize._hungarian import linear_sum_assignment
except ImportError:
    pass

from .common import Benchmark


class Hungarian(Benchmark):
    """
    Test the speed of the Munkres algorithm for solving the
    assignment problem
    """
    params = [
        [100, 300, 500, 700]
    ]

    param_names = ['matrix_size']

    def setup(self, matrix_size):
        nprand.seed(161718)
        self.cost_matrix = rand(matrix_size, matrix_size)

    def time_solve(self, matrix_size):
        linear_sum_assignment(self.cost_matrix)

