"""
benchmarks for the global optimization algorithms

"""
import numpy as np

import time
import inspect
from scipy.optimize import basinhopping, differential_evolution, OptimizeResult

import go_benchmark_functions as gbf
from numpy.testing import Tester, TestCase, assert_array_almost_equal
from collections import defaultdict

NUMTRIALS = 100


class _BenchOptimizers(object):
    """
    A framework for benchmarking optimizers

    Parameters
    ----------
    function_name : string
    function : test object
    minimizer_kwargs : kwargs
        additional keywords passed to the minimizer.  e.g. tol, maxiter
    """

    def __init__(self, function_name, function, **minimizer_kwargs):
        self.function_name = function_name

        self.function = function
        self.fun = function.fun
        if hasattr(function, 'der'):
            self.der = function.der

        self.bounds = function.bounds

        self.minimizer_kwargs = minimizer_kwargs
        self.results = []

    def reset(self):
        self.results = []

    def energy_gradient(self, x):
        return self.fun(x), self.function.der(x)

    # for basinhopping
    def accept_test(self, x_new=None, *args, **kwargs):
        """
        Does the new candidate vector lie inbetween the bounds?

        Returns
        -------
        accept_test : bool
            The candidate vector lies inbetween the bounds
        """
        if not hasattr(self.function, "xmin"):
            return True
        if np.any(x_new < self.function.xmin):
            return False
        if np.any(x_new > self.function.xmax):
            return False
        return True

    def add_result(self, result, t, name):
        """
        Add a result to the list

        Parameters
        ----------
        result : OptimizeResult
            The result from the minimization
        t : float
            Time taken to perform the minimisation
        name : str
            The name of the minimizer used
        """

        result.time = t
        result.name = name
        if not hasattr(result, "njev"):
            result.njev = 0
        if not hasattr(result, "nhev"):
            result.nhev = 0
        self.results.append(result)

    def print_results(self):
        """
        Print the current list of results
        """
        results = self.average_results()
        results = sorted(results, key=lambda x: (x.nfail, x.mean_time))
        print("")
        print("=========================================================")
        print(("Optimizer benchmark: %s" % (self.function_name)))
        print(("dimensions: %d, extra kwargs: %s" %
              (results[0].ndim, str(self.minimizer_kwargs))))
        print(("averaged over %d starting configurations" %
              (results[0].ntrials)))
        print("  Optimizer    nfail   nfev    njev    nhev    time")
        print("---------------------------------------------------------")
        for res in results:
            print(("%30s   | %4d  | %8d  | %4d  | %4d  | %.6g" %
                  (res.name, res.nfail, res.mean_nfev, res.mean_njev,
                   res.mean_nhev, res.mean_time)))

    def average_results(self):
        """
        group the results by minimizer and average over the runs
        """
        grouped_results = defaultdict(list)
        for res in self.results:
            grouped_results[res.name].append(res)

        averaged_results = dict()
        for name, result_list in grouped_results.items():
            newres = OptimizeResult()
            newres.name = name
            newres.mean_nfev = np.mean([r.nfev for r in result_list])
            newres.mean_njev = np.mean([r.njev for r in result_list])
            newres.mean_nhev = np.mean([r.nhev for r in result_list])
            newres.mean_time = np.mean([r.time for r in result_list])
            newres.ntrials = len(result_list)
            newres.nfail = len([r for r in result_list if not r.success])
            try:
                newres.ndim = len(result_list[0].x)
            except TypeError:
                newres.ndim = 1
            averaged_results[name] = newres
        return averaged_results.values()

    def bench_run(self, **minimizer_kwargs):
        """
        Do an optimization test starting at x0 for all the optimizers
        """
        kwargs = self.minimizer_kwargs

        if hasattr(self.fun, "temperature"):
            kwargs["T"] = self.function.temperature
        if hasattr(self.fun, "stepsize"):
            kwargs["stepsize"] = self.function.stepsize
        minimizer_kwargs = {"method": "L-BFGS-B"}
        x0 = self.function.initial_vector()

        # basinhopping - no gradient
        minimizer_kwargs['jac'] = False
        self.function.nfev = 0

        t0 = time.time()

        res = basinhopping(
            self.fun, x0, accept_test=self.accept_test,
            minimizer_kwargs=minimizer_kwargs,
            **kwargs)

        t1 = time.time()
        res.success = self.function.success(res.x)
        res.nfev = self.function.nfev
        self.add_result(res, t1 - t0, 'basinhopping')

        # differential_evolution
        self.function.nfev = 0

        t0 = time.time()

        res = differential_evolution(self.fun,
                                     self.bounds,
                                     popsize=20)

        t1 = time.time()
        res.success = self.function.success(res.x)
        res.nfev = self.function.nfev
        self.add_result(res, t1 - t0, 'differential_evolution')


class BenchGlobalOptimizers(TestCase):
    """
    Benchmark the global optimizers using the go_benchmark_functions
    suite
    """

    def bench_all(self):
        bench_members = inspect.getmembers(gbf, inspect.isclass)
        functions = [item for item in bench_members if
                                    issubclass(item[1], gbf.Benchmark)]

        for name, klass in functions:
            f = klass()
            if name == 'Benchmark':
                continue
            if name is 'LennardJones':
                continue
            if name.startswith('Problem'):
                continue

            b = _BenchOptimizers(name, f)
            for i in range(NUMTRIALS):
                b.bench_run()

            b.print_results()


if __name__ == "__main__":
    Tester().bench(extra_argv=dict())
