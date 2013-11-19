import time
from collections import defaultdict

import numpy as np

import scipy.optimize
from scipy.optimize.optimize import rosen, rosen_der, rosen_hess
import test_functions as funcs


class BenchOptimizers(object):
    """a framework for benchmarking the optimizer
    
    Parameters
    ----------
    function_name : string
    fun : callable
    der : callable
        function that returns the derivitive (jacobian, gradient) of fun
    hess : callable
        function that returns the hessian of fun
    """
    def __init__(self, function_name, fun, der=None, hess=None):
        self.function_name = function_name
        self.fun = fun
        self.der = der
        self.hess = hess
        self.results = []

    def reset(self):
        self.results = []

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
        print("")
        print("---------------------------------------------------------")
        print("Optimizer benchmark: %s" % (self.function_name))
        print("averaged over %d starting configurations" % (results[0].ntrials))
        print("      Optimizer   nfail  nfev   njev   nhev   time")
        print("---------------------------------------------------------")
        for res in results:
            print("%15s   %4d   %4d   %4d   %4d   %.6g" % 
                  (res.name, res.nfail, res.mean_nfev, res.mean_njev, res.mean_nhev, res.mean_time))
    
    def average_results(self):
        """group the results by minimizer and average over the runs"""
        grouped_results = defaultdict(list)
        for res in self.results:
            grouped_results[res.name].append(res)
        
        averaged_results = dict()
        for name, result_list in grouped_results.iteritems():
            newres = scipy.optimize.Result()
            newres.name = name
            newres.mean_nfev = np.mean([r.nfev for r in result_list])
            newres.mean_njev = np.mean([r.njev for r in result_list])
            newres.mean_nhev = np.mean([r.nhev for r in result_list])
            newres.mean_time = np.mean([r.time for r in result_list])
            newres.ntrials = len(result_list)
            newres.nfail = len([r for r in result_list if not r.success])
            averaged_results[name] = newres
        return averaged_results.values()
    
    def bench_run(self, x0):
        """do an optimization test starting at x0 for all the optimizers"""
        kwargs = dict(tol=1e-4)
        
        fonly_methods = ["COBYLA", 'Powell']
        for method in fonly_methods:
            t0 = time.time()
            res = scipy.optimize.minimize(self.fun, x0, method=method, 
                                          **kwargs)
            t1 = time.time()
            self.add_result(res, t1-t0, method)
        
        
        gradient_methods = ['L-BFGS-B', 'BFGS', 'CG', 'TNC', 'SLSQP']
        if self.der is not None:
            for method in gradient_methods:
                t0 = time.time()
                res = scipy.optimize.minimize(self.fun, x0, method=method, 
                                              jac=self.der, **kwargs)
                t1 = time.time()
                self.add_result(res, t1-t0, method)

        hessian_methods = ["Newton-CG", 'dogleg', 'trust-ncg']
        if self.hess is not None:
            for method in hessian_methods:
                t0 = time.time()
                res = scipy.optimize.minimize(self.fun, x0, method=method, 
                                              jac=self.der, hess=self.hess, 
                                              **kwargs)
                t1 = time.time()
                self.add_result(res, t1-t0, method)



def bench_rosenbrock():
    b = BenchOptimizers("Rosenbrock function",
                        fun=rosen, der=rosen_der, hess=rosen_hess)
    
#    # do a single test
#    b.bench_run([0.8, 1.2, 0.7])
#    b.print_results()
    
    # average over multiple starting points
    b.reset()
    for i in xrange(10):
        b.bench_run(np.random.uniform(-3,3,3))
    b.print_results()

def bench_simple_quadratic():
    s = funcs.SimpleQuadratic()
#    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
    b = BenchOptimizers("simple quadratic function",
                        fun=s.fun, der=s.der, hess=s.hess)
    for i in xrange(10):
        b.bench_run(np.random.uniform(-2,2,3))
    b.print_results()

def bench_asymetric_quadratic():
    s = funcs.AsymmetricQuadratic()
#    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
    b = BenchOptimizers("function sum(x**2) + x[0]",
                        fun=s.fun, der=s.der, hess=s.hess)
    for i in xrange(10):
        b.bench_run(np.random.uniform(-2,2,3))
    b.print_results()

def bench_sin_1d():
    fun = lambda x: np.sin(x[0])
    der = lambda x: np.array([np.cos(x[0])])
    b = BenchOptimizers("1d sin function",
                        fun=fun, der=der, hess=None)
    for i in xrange(10):
        b.bench_run(np.random.uniform(-2,2,1))
    b.print_results()

def bench_booth():
    s = funcs.Booth()
#    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
    b = BenchOptimizers("Booth's function",
                        fun=s.fun, der=s.der, hess=None)
    for i in xrange(10):
        b.bench_run(np.random.uniform(0,10,2))
    b.print_results()

def bench_beale():
    s = funcs.Beale()
    print "checking gradient", scipy.optimize.check_grad(s.fun, s.der, np.array([1.1, -2.3]))
    b = BenchOptimizers("Beale's function",
                        fun=s.fun, der=s.der, hess=None)
    for i in xrange(10):
        b.bench_run(np.random.uniform(0,10,2))
    b.print_results()

def bench_LJ():
    s = funcs.LJ()
    print "checking gradient", scipy.optimize.check_grad(s.get_energy, s.get_gradient, np.random.uniform(-2,2,3*4))
    natoms = 4
    b = BenchOptimizers("%d atom Lennard Jones potential" % (natoms),
                        fun=s.get_energy, der=s.get_gradient, hess=None)
    for i in xrange(10):
        b.bench_run(np.random.uniform(-2,2,natoms*3))
    b.print_results()


def main():
    bench_rosenbrock()
    bench_simple_quadratic()
    bench_asymetric_quadratic()
    bench_sin_1d()
    bench_booth()
    bench_beale()
    bench_LJ()

if __name__ == "__main__":
    main()
