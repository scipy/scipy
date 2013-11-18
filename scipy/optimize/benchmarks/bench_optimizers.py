import time

import scipy.optimize
from scipy.optimize.optimize import rosen, rosen_der, rosen_hess

class BenchOptimizers(object):
    def __init__(self):
        self.results = []

    def add_result(self, result, t, name):
        result.time = t
        result.name = name
        if not hasattr(result, "njev"):
            result.njev = 0
        if not hasattr(result, "nhev"):
            result.nhev = 0
        self.results.append(result)
    
    def print_results(self):
        results = sorted(self.results, key=lambda x:x.time)
        print("Optimizer benchmark on the Rosenbrock function sorted by time")
        for res in results:
            if res.success:
                success = "pass"
            else:
                success = "fail"
            print("%15s %s nfev %4d njev %4d nhev %4d time %.6g" % 
                  (res.name, success, res.nfev, res.njev, res.nhev, res.time))
    
    def bench_run(self, fun=rosen, der=rosen_der, hess=rosen_hess, 
                  x0=[0.8, 1.2, 0.7]):
        kwargs = dict(tol=1e-4)
        
        fonly_methods = ["COBYLA", 'Powell']
        for method in fonly_methods:
            t0 = time.time()
            res = scipy.optimize.minimize(fun, x0, method=method, **kwargs)
            t1 = time.time()
            self.add_result(res, t1-t0, method)
        
        
        gradient_methods = ['L-BFGS-B', 'BFGS', 'CG', 'TNC', 'SLSQP']
        for method in gradient_methods:
            t0 = time.time()
            res = scipy.optimize.minimize(fun, x0, method=method, jac=der, 
                                          **kwargs)
            t1 = time.time()
            self.add_result(res, t1-t0, method)

        hessian_methods = ["Newton-CG", 'dogleg', 'trust-ncg']
        for method in hessian_methods:
            t0 = time.time()
            res = scipy.optimize.minimize(fun, x0, method=method, jac=der, 
                                          hess=hess, **kwargs)
            t1 = time.time()
            self.add_result(res, t1-t0, method)


def main():
    b = BenchOptimizers()
    b.bench_run()
    b.print_results()

if __name__ == "__main__":
    main()
