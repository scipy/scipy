from numpy.random import rand

try:
    from scipy.optimize._hungarian import linear_sum_assignment
except ImportError:
    pass

from .common import Benchmark


class Hungarian(Benchmark):
    params = [
        [100, 300, 500, 700]
    ]

    param_names = ['matrix_size']

    def setup(self, matrix_size):
        self.cost_matrix_1 = rand(matrix_size, matrix_size)
        self.cost_matrix_2 = rand(matrix_size, matrix_size)

    def time_solve_minimize(self, matrix_size):
        linear_sum_assignment(self.cost_matrix_1)

    def time_solve_maximize(self, matrix_size):
        linear_sum_assignment(self.cost_matrix_2, maximize=True)

