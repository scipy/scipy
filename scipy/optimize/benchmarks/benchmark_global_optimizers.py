"""
benchmarks for the global optimization algorithms

"""
import numpy as np

import time
from scipy.optimize import basinhopping, differential_evolution, OptimizeResult

import test_functions as funcs
from numpy.testing import Tester, TestCase, assert_array_almost_equal
from collections import defaultdict

NUMTRIALS = 10


class _BenchOptimizers(object):
    """a framework for benchmarking the optimizer

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

        upper_limit = function.xmax
        lower_limit = function.xmin
        limits = np.c_[lower_limit, upper_limit].T
        self.bounds = [(limits[0, idx], limits[1, idx])
                       for idx in range(np.size(limits, 1))]

        self.minimizer_kwargs = minimizer_kwargs
        self.tol = 1e-4

        self.results = []

    def reset(self):
        self.results = []

    def energy_gradient(self, x):
        return self.fun(x), self.function.der(x)

    # for basinhopping
    def get_random_configuration(self):
        if hasattr(self.function, "get_random_configuration"):
            return self.function.get_random_configuration()
        xmin, xmax = self.function.xmin, self.function.xmax
        x = np.random.uniform(xmin[0] + .01, xmax[0] - .01)
        y = np.random.uniform(xmin[1] + .01, xmax[1] - .01)
        return np.array([x, y])

    # for basinhopping
    def accept_test(self, x_new=None, *args, **kwargs):
        if not hasattr(self.function, "xmin"):
            return True
        if np.any(x_new < self.function.xmin):
            return False
        if np.any(x_new > self.function.xmax):
            return False
        return True

    # for basinhopping
    def stop_criterion(self, coords, E, accepted):
        if accepted and E < self.function.target_E + self.tol:
            return True
        else:
            return False

    def found_target(self, res):
        try:
            assert_array_almost_equal(np.abs(res.x),
                                          np.abs(self.function.solution),
                                          decimal=4)
            return True
        except (AssertionError):
            return False

    def add_result(self, result, t, name):
        """
        add a result to the list
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
        print the current list of results
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
        do an optimization test starting at x0 for all the optimizers
        """
        kwargs = self.minimizer_kwargs

        if hasattr(self.fun, "temperature"):
            kwargs["T"] = self.function.temperature
        if hasattr(self.fun, "stepsize"):
            kwargs["stepsize"] = self.function.stepsize
        minimizer_kwargs = {"method": "L-BFGS-B"}
        x0 = self.get_random_configuration()

        # basinhopping - with gradient
        if hasattr(self.function, 'der'):
            minimizer_kwargs['jac'] = True
            t0 = time.time()
            res = basinhopping(
                self.energy_gradient, x0, accept_test=self.accept_test,
                callback=self.stop_criterion, niter=1000,
                minimizer_kwargs=minimizer_kwargs,
                **kwargs)
            t1 = time.time()
            res.success = True
            if not self.found_target(res):
                res.success = False
            self.add_result(res, t1 - t0, 'basinhopping')

        # basinhopping - no gradient
        x0 = self.get_random_configuration()
        minimizer_kwargs['jac'] = False
        t0 = time.time()

        res = basinhopping(
            self.fun, x0, accept_test=self.accept_test,
            callback=self.stop_criterion, niter=1000,
            minimizer_kwargs=minimizer_kwargs,
            **kwargs)

        t1 = time.time()
        res.success = True
        if not self.found_target(res):
            res.success = False
        self.add_result(res, t1 - t0, 'basinhopping - no gradient')

        # differential_evolution
        t0 = time.time()

        res = differential_evolution(self.fun,
                                     self.bounds,
                                     popsize=20,
                                     polish=True)

        t1 = time.time()
        if not self.found_target(res):
            res.success = False
        self.add_result(res, t1 - t0, 'differential_evolution')


class BenchGlobalOptimizers(TestCase):
    """Benchmark the global optimizers"""

    def bench_all(self):
        functions = ['Ackley', 'Levi', 'HolderTable', 'EggHolder',
                     'Schaffer2', 'Schaffer4', 'CrossInTray',
                     'Booth', 'Beale']

        for function in functions:
            f = getattr(funcs, function)()

            b = _BenchOptimizers(function, f)
            for i in range(NUMTRIALS):
                b.bench_run()

            b.print_results()


if __name__ == "__main__":
    Tester().bench(extra_argv=dict())
