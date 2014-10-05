"""
Unit tests for the differential global minimization algorithm.
"""
import numpy as np
import go_benchmark_functions as gbf
import inspect
from numpy.testing import (assert_equal, TestCase, assert_allclose,
                           run_module_suite, assert_almost_equal,
                           assert_string_equal)


class TestDifferentialEvolutionSolver(TestCase):

    def setUp(self):
        self.old_seterr = np.seterr(invalid='raise')
        bench_members = inspect.getmembers(gbf, inspect.isclass)

        self.benchmark_functions = [item for item in bench_members if
                                    issubclass(item[1], gbf.Benchmark)]

    def tearDown(self):
        np.seterr(**self.old_seterr)


    def test_optimum_solution(self):
        """
        Check that the function returns the global minimum if given
        the optimal solution
        """
        for name, klass in self.benchmark_functions:
            f = klass()
            assert(f.success(f.global_optimum))


if __name__ == '__main__':
    run_module_suite()
