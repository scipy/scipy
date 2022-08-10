"""
Unit tests for the global optimization benchmark functions
"""
import numpy as np
from .. import go_benchmark_functions as gbf
import inspect


class TestGoBenchmarkFunctions:

    def setup_method(self):
        bench_members = inspect.getmembers(gbf, inspect.isclass)
        self.benchmark_functions = {it[0]:it[1] for it in bench_members if
                                    issubclass(it[1], gbf.Benchmark)}

    def teardown_method(self):
        pass

    def test_optimum_solution(self):
        # Check that the function returns the global minimum if given
        # the optimal solution
        for name, klass in self.benchmark_functions.items():
            # LennardJones is filtered here because there are many global
            # optimima that give the same minimum energy
            if (name in ['Benchmark', 'LennardJones'] or
                 name.startswith('Problem')):
                continue

            f = klass()

            if name in ['Damavandi', 'Csendes']:
                with np.errstate(divide='ignore', invalid='ignore'):
                    print(name, f.fun(np.asarray(f.global_optimum[0])),
                          f.fglob)
                    assert np.isnan(f.fun(np.asarray(f.global_optimum[0])))
                    continue

            print(name, f.fun(np.asarray(f.global_optimum[0])), f.fglob)
            assert f.success(f.global_optimum[0])

    def test_solution_exists(self):
        # Every benchmark function should have a minimum energy
        for name, klass in self.benchmark_functions.items():
            if name == 'Benchmark':
                continue

            f = klass()
            # should result in an attribute error if it doesn't exist
            val = f.fglob

    def test_bounds_access_subscriptable(self):
        # In Python 2 zip returns a list which is subscriptable
        # In Python 3 zip returns a zip object, which is not subscriptable
        for name, klass in self.benchmark_functions.items():
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
        assert len(x0) == 10

        # the bounds should be the correct length now.
        bounds = L.bounds
        assert len(bounds) == 10

        assert L.N == 10
