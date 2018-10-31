from __future__ import division, print_function, absolute_import

import os
import time
import inspect
import json
import traceback
from collections import defaultdict, OrderedDict

import numpy as np

try:
    import scipy.optimize
    from scipy.optimize.optimize import rosen, rosen_der, rosen_hess
    from scipy.optimize import (leastsq, basinhopping, differential_evolution,
                                OptimizeResult)
except ImportError:
    pass

from . import test_functions as funcs
from . import go_benchmark_functions as gbf
from .common import Benchmark
from .lsq_problems import extract_lsq_problems


class _BenchOptimizers(Benchmark):
    """a framework for benchmarking the optimizer

    Parameters
    ----------
    function_name : string
    fun : callable
    der : callable
        function that returns the derivative (jacobian, gradient) of fun
    hess : callable
        function that returns the hessian of fun
    minimizer_kwargs : kwargs
        additional keywords passed to the minimizer.  e.g. tol, maxiter
    """
    def __init__(self, function_name, fun, der=None, hess=None,
                 **minimizer_kwargs):
        self.function_name = function_name
        self.fun = fun
        self.der = der
        self.hess = hess
        self.minimizer_kwargs = minimizer_kwargs
        if "tol" not in minimizer_kwargs:
            minimizer_kwargs["tol"] = 1e-4

        self.results = []

    @classmethod
    def from_funcobj(cls, function_name, function, **minimizer_kwargs):
        self = cls.__new__(cls)
        self.function_name = function_name

        self.function = function
        self.fun = function.fun
        if hasattr(function, 'der'):
            self.der = function.der

        self.bounds = function.bounds

        self.minimizer_kwargs = minimizer_kwargs
        self.results = []
        return self

    def reset(self):
        self.results = []

    def energy_gradient(self, x):
        return self.fun(x), self.function.der(x)

    def add_result(self, result, t, name):
        """add a result to the list"""
        result.time = t
        result.name = name
        if not hasattr(result, "njev"):
            result.njev = 0
        if not hasattr(result, "nhev"):
            result.nhev = 0
        self.results.append(result)

    def print_results(self):
        """print the current list of results"""
        results = self.average_results()
        results = sorted(results, key=lambda x: (x.nfail, x.mean_time))
        if not results:
            return
        print("")
        print("=========================================================")
        print("Optimizer benchmark: %s" % (self.function_name))
        print("dimensions: %d, extra kwargs: %s" % (results[0].ndim, str(self.minimizer_kwargs)))
        print("averaged over %d starting configurations" % (results[0].ntrials))
        print("  Optimizer    nfail   nfev    njev    nhev    time")
        print("---------------------------------------------------------")
        for res in results:
            print("%11s  | %4d  | %4d  | %4d  | %4d  | %.6g" %
                  (res.name, res.nfail, res.mean_nfev, res.mean_njev, res.mean_nhev, res.mean_time))

    def average_results(self):
        """group the results by minimizer and average over the runs"""
        grouped_results = defaultdict(list)
        for res in self.results:
            grouped_results[res.name].append(res)

        averaged_results = dict()
        for name, result_list in grouped_results.items():
            newres = scipy.optimize.OptimizeResult()
            newres.name = name
            newres.mean_nfev = np.mean([r.nfev for r in result_list])
            newres.mean_njev = np.mean([r.njev for r in result_list])
            newres.mean_nhev = np.mean([r.nhev for r in result_list])
            newres.mean_time = np.mean([r.time for r in result_list])
            newres.ntrials = len(result_list)
            newres.nfail = len([r for r in result_list if not r.success])
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
        self.add_result(res, t1 - t0, 'basinh.')

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
        self.add_result(res, t1 - t0, 'DE')

    def bench_run_global(self, numtrials=50, methods=None):
        """
        Run the optimization tests for the required minimizers.
        """

        if methods is None:
            methods = ['DE', 'basinh.']

        method_fun = {'DE': self.run_differentialevolution,
                      'basinh.': self.run_basinhopping}

        for i in range(numtrials):
            for m in methods:
                method_fun[m]()

    def bench_run(self, x0, methods=None, **minimizer_kwargs):
        """do an optimization test starting at x0 for all the optimizers"""
        kwargs = self.minimizer_kwargs

        if methods is None:
            methods = ["COBYLA", 'Powell', 'nelder-mead',
                       'L-BFGS-B', 'BFGS', 'CG', 'TNC', 'SLSQP',
                       "Newton-CG", 'dogleg', 'trust-ncg', 'trust-exact',
                       'trust-krylov']

        # L-BFGS-B, BFGS can use gradients, but examine performance when
        # numerical differentiation is used.
        fonly_methods = ["COBYLA", 'Powell', 'nelder-mead', 'L-BFGS-B', 'BFGS']
        for method in fonly_methods:
            if method not in methods:
                continue
            t0 = time.time()
            res = scipy.optimize.minimize(self.fun, x0, method=method,
                                          **kwargs)
            t1 = time.time()
            self.add_result(res, t1-t0, method)

        gradient_methods = ['L-BFGS-B', 'BFGS', 'CG', 'TNC', 'SLSQP']
        if self.der is not None:
            for method in gradient_methods:
                if method not in methods:
                    continue
                t0 = time.time()
                res = scipy.optimize.minimize(self.fun, x0, method=method,
                                              jac=self.der, **kwargs)
                t1 = time.time()
                self.add_result(res, t1-t0, method)

        hessian_methods = ["Newton-CG", 'dogleg', 'trust-ncg',
                           'trust-exact', 'trust-krylov']
        if self.hess is not None:
            for method in hessian_methods:
                if method not in methods:
                    continue
                t0 = time.time()
                res = scipy.optimize.minimize(self.fun, x0, method=method,
                                              jac=self.der, hess=self.hess,
                                              **kwargs)
                t1 = time.time()
                self.add_result(res, t1-t0, method)


class BenchSmoothUnbounded(Benchmark):
    """Benchmark the optimizers with smooth, unbounded, functions"""
    params = [
        ['rosenbrock_slow', 'rosenbrock_nograd', 'rosenbrock', 'rosenbrock_tight',
         'simple_quadratic', 'asymmetric_quadratic',
         'sin_1d', 'booth', 'beale', 'LJ'],
        ["COBYLA", 'Powell', 'nelder-mead',
         'L-BFGS-B', 'BFGS', 'CG', 'TNC', 'SLSQP',
         "Newton-CG", 'dogleg', 'trust-ncg', 'trust-exact',
         'trust-krylov'],
        ["mean_nfev", "mean_time"]
    ]
    param_names = ["test function", "solver", "result type"]

    def setup(self, func_name, method_name, ret_val):
        b = getattr(self, 'run_' + func_name)(methods=[method_name])
        r = b.average_results().get(method_name)
        if r is None:
            raise NotImplementedError()
        self.result = getattr(r, ret_val)

    def track_all(self, func_name, method_name, ret_val):
        return self.result

    # SlowRosen has a 50us delay on each function evaluation. By comparing to
    # rosenbrock_nograd it should be possible to figure out how much time a
    # minimizer uses internally, compared to the time required for function
    # evaluation.
    def run_rosenbrock_slow(self, methods=None):
        s = funcs.SlowRosen()
        b = _BenchOptimizers("Rosenbrock function",
                             fun=s.fun)
        for i in range(10):
            b.bench_run(np.random.uniform(-3, 3, 3), methods=methods)
        return b

    # see what the performance of the solvers are if numerical differentiation
    # has to be used.
    def run_rosenbrock_nograd(self, methods=None):
        b = _BenchOptimizers("Rosenbrock function",
                             fun=rosen)
        for i in range(10):
            b.bench_run(np.random.uniform(-3, 3, 3), methods=methods)
        return b

    def run_rosenbrock(self, methods=None):
        b = _BenchOptimizers("Rosenbrock function",
                             fun=rosen, der=rosen_der, hess=rosen_hess)
        for i in range(10):
            b.bench_run(np.random.uniform(-3, 3, 3), methods=methods)
        return b

    def run_rosenbrock_tight(self, methods=None):
        b = _BenchOptimizers("Rosenbrock function",
                             fun=rosen, der=rosen_der, hess=rosen_hess,
                             tol=1e-8)
        for i in range(10):
            b.bench_run(np.random.uniform(-3, 3, 3), methods=methods)
        return b

    def run_simple_quadratic(self, methods=None):
        s = funcs.SimpleQuadratic()
        #    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("simple quadratic function",
                             fun=s.fun, der=s.der, hess=s.hess)
        for i in range(10):
            b.bench_run(np.random.uniform(-2, 2, 3), methods=methods)
        return b

    def run_asymmetric_quadratic(self, methods=None):
        s = funcs.AsymmetricQuadratic()
        #    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("function sum(x**2) + x[0]",
                             fun=s.fun, der=s.der, hess=s.hess)
        for i in range(10):
            b.bench_run(np.random.uniform(-2, 2, 3), methods=methods)
        return b

    def run_sin_1d(self, methods=None):
        fun = lambda x: np.sin(x[0])
        der = lambda x: np.array([np.cos(x[0])])
        b = _BenchOptimizers("1d sin function",
                             fun=fun, der=der, hess=None)
        for i in range(10):
            b.bench_run(np.random.uniform(-2, 2, 1), methods=methods)
        return b

    def run_booth(self, methods=None):
        s = funcs.Booth()
        #    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("Booth's function",
                             fun=s.fun, der=s.der, hess=None)
        for i in range(10):
            b.bench_run(np.random.uniform(0, 10, 2), methods=methods)
        return b

    def run_beale(self, methods=None):
        s = funcs.Beale()
        #    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("Beale's function",
                             fun=s.fun, der=s.der, hess=None)
        for i in range(10):
            b.bench_run(np.random.uniform(0, 10, 2), methods=methods)
        return b

    def run_LJ(self, methods=None):
        s = funcs.LJ()
        # print "checking gradient", scipy.optimize.check_grad(s.get_energy, s.get_gradient,
        # np.random.uniform(-2,2,3*4))
        natoms = 4
        b = _BenchOptimizers("%d atom Lennard Jones potential" % (natoms),
                             fun=s.fun, der=s.der, hess=None)
        for i in range(10):
            b.bench_run(np.random.uniform(-2, 2, natoms*3), methods=methods)
        return b


class BenchLeastSquares(Benchmark):
    """Class for benchmarking nonlinear least squares solvers."""
    problems = extract_lsq_problems()
    params = [
        list(problems.keys()),
        ["average time", "nfev", "success"]
    ]
    param_names = [
        "problem", "result type"
    ]

    def track_all(self, problem_name, result_type):
        problem = self.problems[problem_name]

        if problem.lb is not None or problem.ub is not None:
            raise NotImplementedError

        ftol = 1e-5

        if result_type == 'average time':
            n_runs = 10
            t0 = time.time()
            for _ in range(n_runs):
                leastsq(problem.fun, problem.x0, Dfun=problem.jac, ftol=ftol,
                        full_output=True)
            return (time.time() - t0) / n_runs

        x, cov_x, info, message, ier = leastsq(
            problem.fun, problem.x0, Dfun=problem.jac,
            ftol=ftol, full_output=True
        )
        if result_type == 'nfev':
            return info['nfev']
        elif result_type == 'success':
            return int(problem.check_answer(x, ftol))
        else:
            raise NotImplementedError


try:
    # the value of SCIPY_XSLOW is used to control how many repeats of each
    # function
    slow = int(os.environ.get('SCIPY_XSLOW', 0))
except ValueError:
    pass


_func_names = os.environ.get('SCIPY_GLOBAL_BENCH', [])
if _func_names:
    if not slow:
        slow = 100
    _func_names = [x.strip() for x in _func_names.split(',')]


class BenchGlobal(Benchmark):
    """
    Benchmark the global optimizers using the go_benchmark_functions
    suite
    """
    timeout = 300

    _functions = OrderedDict([
        item for item in inspect.getmembers(gbf, inspect.isclass)
        if (issubclass(item[1], gbf.Benchmark) and
            item[0] not in ('Benchmark') and
            not item[0].startswith('Problem'))
    ])

    if _func_names:
        _filtered_funcs = OrderedDict()
        for name in _func_names:
            if name in _functions:
                _filtered_funcs[name] = _functions.get(name)
        _functions = _filtered_funcs

    if not slow:
        _functions = {'AMGM': None}

    params = [
        list(_functions.keys()),
        ["success%", "<nfev>"],
        ['DE', 'basinh.'],
    ]
    param_names = ["test function", "result type", "solver"]

    def __init__(self):
        self.enabled = bool(slow)
        self.numtrials = slow

        self.dump_fn = os.path.join(os.path.dirname(__file__), '..', 'global-bench-results.json')
        self.results = {}

    def setup(self, name, ret_value, solver):
        if not self.enabled:
            print("BenchGlobal.track_all not enabled --- export SCIPY_XSLOW=slow to enable,\n"
                  "'slow' iterations of each benchmark will be run.\n"
                  "Note that it can take several hours to run; intermediate output\n"
                  "can be found under benchmarks/global-bench-results.json\n"
                  "You can specify functions to benchmark via SCIPY_GLOBAL_BENCH=AMGM,Adjiman,...")
            raise NotImplementedError()
        # load json backing file
        with open(self.dump_fn, 'r') as f:
            self.results = json.load(f)

    def teardown(self, name, ret_value, solver):
        with open(self.dump_fn, 'w') as f:
            json.dump(self.results, f, indent=2, sort_keys=True)

    def track_all(self, name, ret_value, solver):
        if name in self.results and solver in self.results[name]:
            # have we done the function, and done the solver?
            # if so, then just return the ret_value
            av_results = self.results[name]
            if ret_value == 'success%':
                return 100 * av_results[solver]['nsuccess'] / av_results[solver]['ntrials']
            elif ret_value == '<nfev>':
                return av_results[solver]['mean_nfev']
            else:
                raise ValueError()

        klass = self._functions[name]
        f = klass()
        try:
            b = _BenchOptimizers.from_funcobj(name, f)
            with np.errstate(all='ignore'):
                b.bench_run_global(methods=[solver],
                                   numtrials=self.numtrials)

            av_results = b.average_results()

            if name not in self.results:
                self.results[name] = {}
            self.results[name][solver] = av_results[solver]

            if ret_value == 'success%':
                return 100 * av_results[solver]['nsuccess'] / av_results[solver]['ntrials']
            elif ret_value == '<nfev>':
                return av_results[solver]['mean_nfev']
            else:
                raise ValueError()
        except Exception:
            print("".join(traceback.format_exc()))
            self.results[name] = "".join(traceback.format_exc())

    def setup_cache(self):
        if not self.enabled:
            return

        # create the logfile to start with
        with open(self.dump_fn, 'w') as f:
            json.dump({}, f, indent=2)
