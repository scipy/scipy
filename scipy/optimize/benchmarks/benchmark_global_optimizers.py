from __future__ import division
"""
benchmarks for the global optimization algorithms

"""
import numpy as np

import time
import inspect
from scipy.optimize import basinhopping, differential_evolution, OptimizeResult
from scipy.optimize import anneal

import go_benchmark_functions as gbf
from numpy.testing import *
from collections import defaultdict

NUMTRIALS = 50


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
            newres.nsuccess = len([r for r in result_list if r.success])
            try:
                newres.ndim = len(result_list[0].x)
            except TypeError:
                newres.ndim = 1
            averaged_results[name] = newres
        return averaged_results

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

    def run_basinhopping(self):
        """
        Do an optimization run for basinhopping
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

    def run_differentialevolution(self):
        """
        Do an optimization run for differential_evolution
        """
        self.function.nfev = 0

        t0 = time.time()

        res = differential_evolution(self.fun,
                                     self.bounds,
                                     popsize=20)

        t1 = time.time()
        res.success = self.function.success(res.x)
        res.nfev = self.function.nfev
        self.add_result(res, t1 - t0, 'differential_evolution')

    def run_anneal(self):
        """
        Do an optimization run for simulated annealing
        """
        self.function.nfev = 0
        x0 = self.function.initial_vector()

        t0 = time.time()

        result = anneal(self.fun, x0, disp=False, lower=self.function.xmin,
                        upper=self.function.xmax)

        t1 = time.time()

        res = OptimizeResult()
        res.x = result[0]
        res.success = self.function.success(res.x)
        res.nfev = self.function.nfev
        self.add_result(res, t1 - t0, 'anneal')

    def bench_run(self, trials=NUMTRIALS):
        """
        Run the optimization tests for the required minimizers.
        """
        for i in range(trials):
            self.run_differentialevolution()
            self.run_basinhopping()
            self.run_anneal()


class BenchGlobalOptimizers(TestCase):
    """
    Benchmark the global optimizers using the go_benchmark_functions
    suite
    """

    def bench_all(self):
        bench_members = inspect.getmembers(gbf, inspect.isclass)
        functions = [item for item in bench_members if
                                    issubclass(item[1], gbf.Benchmark)]

        print('')
        print('------------------------------')
        print('Benchmarking Global Optimizers')
        print('------------------------------')
        print('Trials: {0}'.format(NUMTRIALS))
        print('{0:20} {1:^30} {2:^30}'.format('',
                                              'success %',
                                              'nfev'))

        print('{0:20} {1:>10} {2:>10} {3:>10} {4:>10} {5:>10}'
              ' {6:>10}'.format('Problem',
                                'basin %',
                                'diffev %',
                                'anneal %',
                                'basin',
                                'diffev',
                                'anneal'))

        for name, klass in functions:
            if name == 'Benchmark':
                continue
            if name == 'LennardJones':
                continue
            if name.startswith('Problem'):
                continue

            f = klass()
            b = _BenchOptimizers(name, f)
            b.bench_run()
            av_results = b.average_results()

            print('{0:20} {1:>10.1f} {2:>10.1f} {3:>10.1f} {4:>10d} {5:>10d}'
                  ' {6:>10d}'.format(name[0: 20],
                100 * av_results['basinhopping'].nsuccess / NUMTRIALS,
                100 * av_results['differential_evolution'].nsuccess / NUMTRIALS,
                100 * av_results['anneal'].nsuccess / NUMTRIALS,
                int(av_results['basinhopping'].mean_nfev),
                int(av_results['differential_evolution'].mean_nfev),
                int(av_results['anneal'].mean_nfev)))


if __name__ == "__main__":
    Tester().bench(extra_argv=dict())
