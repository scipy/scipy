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
            # LennardJones is filtered here because there are many global
            # optimima that give the same minimum energy
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

    def test_bounds_access_subscriptable(self):
        # In Python 2 zip returns a list which is subscriptable
        # In Python 3 zip returns a zip object, which is not subscriptable
        for name, klass in self.benchmark_functions:
            if (name == 'Benchmark' or name.startswith('Problem')):
                continue

            f = klass()
            bounds = f.bounds[0]

    def test_redimension(self):
        # check that problems can be redimensioned, use LJ for this.
        LJ = self.benchmark_functions['LennardJones']
        L = LJ()
        L.change_dimensions(10)

        # if we change the size of the problem then the initial vector has to
        # resize
        x0 = L.initial_vector()
        assert_(len(x0) == 10)

        # the bounds should be the correct length now.
        bounds = L.bounds
        assert_(len(bounds) == 10)

        assert_(L.N == 10)


if __name__ == '__main__':
    run_module_suite()
