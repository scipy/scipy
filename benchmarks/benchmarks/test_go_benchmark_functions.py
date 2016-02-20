"""
Unit tests for the global optimization benchmark functions
"""
from __future__ import division, print_function, absolute_import
import numpy as np
from . import go_benchmark_functions as gbf
import inspect
from numpy.testing import TestCase, run_module_suite, assert_


class TestGoBenchmarkFunctions(TestCase):

    def setUp(self):
        bench_members = inspect.getmembers(gbf, inspect.isclass)
        self.benchmark_functions = [item for item in bench_members if
                                    issubclass(item[1], gbf.Benchmark)]

    def tearDown(self):
        pass

    def test_optimum_solution(self):
        # Check that the function returns the global minimum if given
        # the optimal solution
        for name, klass in self.benchmark_functions:
            if (name in ['Benchmark', 'LennardJones'] or
                 name.startswith('Problem')):
                continue

            f = klass()
            print(name, f.fun(np.asarray(f.global_optimum[0])), f.fglob)
            assert_(f.success(f.global_optimum[0]))

    def test_solution_exists(self):
        # Every benchmark function should have a minimum energy
        for name, klass in self.benchmark_functions:
            if name == 'Benchmark':
                continue

            f = klass()
            # should result in an attribute error if it doesn't exist
            val = f.fglob


if __name__ == '__main__':
    run_module_suite()
