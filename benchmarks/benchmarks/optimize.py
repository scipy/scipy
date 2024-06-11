import os
import time
import inspect
import json
import traceback
from collections import defaultdict

import numpy as np

from . import test_functions as funcs
from . import go_benchmark_functions as gbf
from .common import Benchmark, is_xslow, safe_import
from .lsq_problems import extract_lsq_problems

with safe_import():
    import scipy.optimize
    from scipy.optimize.optimize import rosen, rosen_der, rosen_hess
    from scipy.optimize import (leastsq, basinhopping, differential_evolution,
                                dual_annealing, shgo, direct)
    from scipy.optimize._minimize import MINIMIZE_METHODS
    from .cutest.calfun import calfun
    from .cutest.dfoxs import dfoxs


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
        print("dimensions: %d, extra kwargs: %s" %
              (results[0].ndim, str(self.minimizer_kwargs)))
        print("averaged over %d starting configurations" % (results[0].ntrials))
        print("  Optimizer    nfail   nfev    njev    nhev    time")
        print("---------------------------------------------------------")
        for res in results:
            print("%11s  | %4d  | %4d  | %4d  | %4d  | %.6g" %
                  (res.name, res.nfail, res.mean_nfev,
                   res.mean_njev, res.mean_nhev, res.mean_time))

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
            funs = [r.fun for r in result_list]
            newres.max_obj = np.max(funs)
            newres.min_obj = np.min(funs)
            newres.mean_obj = np.mean(funs)

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
        Does the new candidate vector lie in between the bounds?

        Returns
        -------
        accept_test : bool
            The candidate vector lies in between the bounds
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

    def run_direct(self):
        """
        Do an optimization run for direct
        """
        self.function.nfev = 0

        t0 = time.time()

        res = direct(self.fun,
                     self.bounds)

        t1 = time.time()
        res.success = self.function.success(res.x)
        res.nfev = self.function.nfev
        self.add_result(res, t1 - t0, 'DIRECT')

    def run_shgo(self):
        """
        Do an optimization run for shgo
        """
        self.function.nfev = 0

        t0 = time.time()

        res = shgo(self.fun,
                   self.bounds)

        t1 = time.time()
        res.success = self.function.success(res.x)
        res.nfev = self.function.nfev
        self.add_result(res, t1 - t0, 'SHGO')

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

    def run_dualannealing(self):
        """
        Do an optimization run for dual_annealing
        """
        self.function.nfev = 0

        t0 = time.time()

        res = dual_annealing(self.fun,
                             self.bounds)

        t1 = time.time()
        res.success = self.function.success(res.x)
        res.nfev = self.function.nfev
        self.add_result(res, t1 - t0, 'DA')

    def bench_run_global(self, numtrials=50, methods=None):
        """
        Run the optimization tests for the required minimizers.
        """

        if methods is None:
            methods = ['DE', 'basinh.', 'DA', 'DIRECT', 'SHGO']

        stochastic_methods = ['DE', 'basinh.', 'DA']

        method_fun = {'DE': self.run_differentialevolution,
                      'basinh.': self.run_basinhopping,
                      'DA': self.run_dualannealing,
                      'DIRECT': self.run_direct,
                      'SHGO': self.run_shgo, }

        for m in methods:
            if m in stochastic_methods:
                for i in range(numtrials):
                    method_fun[m]()
            else:
                method_fun[m]()

    def bench_run(self, x0, methods=None, **minimizer_kwargs):
        """do an optimization test starting at x0 for all the optimizers"""
        kwargs = self.minimizer_kwargs

        if methods is None:
            methods = MINIMIZE_METHODS

        # L-BFGS-B, BFGS, trust-constr, SLSQP can use gradients, but examine
        # performance when numerical differentiation is used.
        fonly_methods = ["COBYLA", 'COBYQA', 'Powell', 'nelder-mead',
                         'L-BFGS-B', 'BFGS', 'trust-constr', 'SLSQP']
        for method in fonly_methods:
            if method not in methods:
                continue
            t0 = time.time()
            res = scipy.optimize.minimize(self.fun, x0, method=method,
                                          **kwargs)
            t1 = time.time()
            self.add_result(res, t1-t0, method)

        gradient_methods = ['L-BFGS-B', 'BFGS', 'CG', 'TNC', 'SLSQP',
                            'trust-constr']
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
                           'trust-exact', 'trust-krylov', 'trust-constr']
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
        ["COBYLA", 'COBYQA', 'Powell', 'nelder-mead',
         'L-BFGS-B', 'BFGS', 'CG', 'TNC', 'SLSQP',
         "Newton-CG", 'dogleg', 'trust-ncg', 'trust-exact',
         'trust-krylov', 'trust-constr'],
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
        #    print "checking gradient",
        #    scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("simple quadratic function",
                             fun=s.fun, der=s.der, hess=s.hess)
        for i in range(10):
            b.bench_run(np.random.uniform(-2, 2, 3), methods=methods)
        return b

    def run_asymmetric_quadratic(self, methods=None):
        s = funcs.AsymmetricQuadratic()
        #    print "checking gradient",
        #    scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("function sum(x**2) + x[0]",
                             fun=s.fun, der=s.der, hess=s.hess)
        for i in range(10):
            b.bench_run(np.random.uniform(-2, 2, 3), methods=methods)
        return b

    def run_sin_1d(self, methods=None):
        def fun(x):
            return np.sin(x[0])

        def der(x):
            return np.array([np.cos(x[0])])

        b = _BenchOptimizers("1d sin function",
                             fun=fun, der=der, hess=None)
        for i in range(10):
            b.bench_run(np.random.uniform(-2, 2, 1), methods=methods)
        return b

    def run_booth(self, methods=None):
        s = funcs.Booth()
        #    print "checking gradient",
        #    scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("Booth's function",
                             fun=s.fun, der=s.der, hess=None)
        for i in range(10):
            b.bench_run(np.random.uniform(0, 10, 2), methods=methods)
        return b

    def run_beale(self, methods=None):
        s = funcs.Beale()
        #    print "checking gradient",
        #    scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
        b = _BenchOptimizers("Beale's function",
                             fun=s.fun, der=s.der, hess=None)
        for i in range(10):
            b.bench_run(np.random.uniform(0, 10, 2), methods=methods)
        return b

    def run_LJ(self, methods=None):
        s = funcs.LJ()
        # print "checking gradient",
        # scipy.optimize.check_grad(s.get_energy, s.get_gradient,
        #                           np.random.uniform(-2,2,3*4))
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


# `export SCIPY_XSLOW=1` to enable BenchGlobal.track_all
# `export SCIPY_GLOBAL_BENCH=AMGM,Adjiman,...` to run specific tests
# `export SCIPY_GLOBAL_BENCH_NUMTRIALS=10` to specify n_iterations, default 100
#
# Note that it can take several hours to run; intermediate output
# can be found under benchmarks/global-bench-results.json


class BenchGlobal(Benchmark):
    """
    Benchmark the global optimizers using the go_benchmark_functions
    suite
    """
    timeout = 300

    _functions = dict([
        item for item in inspect.getmembers(gbf, inspect.isclass)
        if (issubclass(item[1], gbf.Benchmark) and
            item[0] not in ('Benchmark') and
            not item[0].startswith('Problem'))
    ])

    if not is_xslow():
        _enabled_functions = []
    elif 'SCIPY_GLOBAL_BENCH' in os.environ:
        _enabled_functions = [x.strip() for x in
                              os.environ['SCIPY_GLOBAL_BENCH'].split(',')]
    else:
        _enabled_functions = list(_functions.keys())

    params = [
        list(_functions.keys()),
        ["success%", "<nfev>", "average time"],
        ['DE', 'basinh.', 'DA', 'DIRECT', 'SHGO'],
    ]
    param_names = ["test function", "result type", "solver"]

    def __init__(self):
        self.enabled = is_xslow()
        try:
            self.numtrials = int(os.environ['SCIPY_GLOBAL_BENCH_NUMTRIALS'])
        except (KeyError, ValueError):
            self.numtrials = 100

        self.dump_fn = os.path.join(os.path.dirname(__file__),
                                    '..',
                                    'global-bench-results.json',)
        self.results = {}

    def setup(self, name, ret_value, solver):
        if name not in self._enabled_functions:
            raise NotImplementedError("skipped")

        # load json backing file
        with open(self.dump_fn) as f:
            self.results = json.load(f)

    def teardown(self, name, ret_value, solver):
        if not self.enabled:
            return

        with open(self.dump_fn, 'w') as f:
            json.dump(self.results, f, indent=2, sort_keys=True)

    def track_all(self, name, ret_value, solver):
        if name in self.results and solver in self.results[name]:
            # have we done the function, and done the solver?
            # if so, then just return the ret_value
            av_results = self.results[name]
            if ret_value == 'success%':
                return (100 * av_results[solver]['nsuccess']
                        / av_results[solver]['ntrials'])
            elif ret_value == '<nfev>':
                return av_results[solver]['mean_nfev']
            elif ret_value == 'average time':
                return av_results[solver]['mean_time']
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
                return (100 * av_results[solver]['nsuccess']
                        / av_results[solver]['ntrials'])
            elif ret_value == '<nfev>':
                return av_results[solver]['mean_nfev']
            elif ret_value == 'average time':
                return av_results[solver]['mean_time']
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


class BenchDFO(Benchmark):
    """
    Benchmark the optimizers with the CUTEST DFO benchmark of Moré and Wild.
    The original benchmark suite is available at
    https://github.com/POptUS/BenDFO
    """

    params = [
        list(range(53)),  # adjust which problems to solve
        ["COBYLA", "COBYQA", "SLSQP", "Powell", "nelder-mead", "L-BFGS-B",
         "BFGS",
         "trust-constr"],  # note: methods must also be listed in bench_run
        ["mean_nfev", "min_obj"],  # defined in average_results
    ]
    param_names = ["DFO benchmark problem number", "solver", "result type"]

    def setup(self, prob_number, method_name, ret_val):
        probs = np.loadtxt(os.path.join(os.path.dirname(__file__),
                                        "cutest", "dfo.txt"))
        params = probs[prob_number]
        nprob = int(params[0])
        n = int(params[1])
        m = int(params[2])
        s = params[3]
        factor = 10 ** s

        def func(x):
            return calfun(x, m, nprob)

        x0 = dfoxs(n, nprob, factor)
        b = getattr(self, "run_cutest")(
            func, x0, prob_number=prob_number, methods=[method_name]
        )
        r = b.average_results().get(method_name)
        if r is None:
            raise NotImplementedError()
        self.result = getattr(r, ret_val)

    def track_all(self, prob_number, method_name, ret_val):
        return self.result

    def run_cutest(self, func, x0, prob_number, methods=None):
        if methods is None:
            methods = MINIMIZE_METHODS
        b = _BenchOptimizers(f"DFO benchmark problem {prob_number}", fun=func)
        b.bench_run(x0, methods=methods)
        return b
