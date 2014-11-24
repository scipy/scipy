"""
Unit tests for the global optimization benchmark functions
"""
import numpy as np
import go_benchmark_functions as gbf
import inspect
from numpy.testing import (assert_equal, TestCase, assert_allclose,
                           run_module_suite, assert_almost_equal,
                           assert_string_equal)


class TestGoBenchmarkFunctions(TestCase):

    def setUp(self):
        #self.old_seterr = np.seterr(invalid='raise')

        bench_members = inspect.getmembers(gbf, inspect.isclass)
        self.benchmark_functions = [item for item in bench_members if
                                    issubclass(item[1], gbf.Benchmark)]

    def tearDown(self):
        #np.seterr(**self.old_seterr)
        pass

    def test_optimum_solution(self):
        #Check that the function returns the global minimum if given
        #the optimal solution
        for name, klass in self.benchmark_functions:
            if name == 'Benchmark':
                continue
            if name == 'LennardJones':
                continue
            if name.startswith('Problem'):
                continue

            f = klass()
            print name, f.fun(np.asarray(f.global_optimum[0])), f.fglob
            assert(f.success(f.global_optimum[0]))

    def test_solution_exists(self):
        #Every benchmark function should have a minimum energy
        for name, klass in self.benchmark_functions:
            if name == 'Benchmark':
                continue

            f = klass()
            #should result in an attribute error if it doesn't exist
            val = f.fglob


if __name__ == '__main__':
    run_module_suite()
