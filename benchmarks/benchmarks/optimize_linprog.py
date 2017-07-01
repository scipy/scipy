"""
Benchmarks for Linear Programming
"""
from __future__ import division, print_function, absolute_import

# Import testing parameters
try:
    from scipy.optimize import linprog
    from scipy.linalg import toeplitz
    from scipy.optimize.tests.test_linprog import lpgen_2d
    import numpy as np
except ImportError:
    pass

from .common import Benchmark

methods = ["simplex","interior-point"]
def klee_minty(D):
    A_1 = np.array([2**(i+1) if i > 0 else 1 for i in range(D)])
    A1_ = np.zeros(D)
    A1_[0] = 1
    A_ub = toeplitz(A_1,A1_)
    b_ub = np.array([5**(i+1) for i in range(D)])
    c = -np.array([2**(D-i-1) for i in range(D)])
    return c, A_ub, b_ub

class KleeMinty(Benchmark):

    params = [
        methods,
        [3,6,9]
    ]
    param_names = ['method', 'dimensions']

    def setup(self, meth, dims):
        self.c, self.A_ub, self.b_ub = klee_minty(dims)
        self.meth = meth

    def time_klee_minty(self, meth, dims):
        linprog(c=self.c, A_ub=self.A_ub, b_ub=self.b_ub, 
                method=self.meth)
        
class LpGen(Benchmark):
    params = [
        methods,
        range(20,100,20),
        range(20,100,20)        
    ]
    param_names = ['method','m','n']
    
    def setup(self, meth, m, n):
        self.A,self.b,self.c = lpgen_2d(m,n)
        self.meth = meth
    
    def time_lpgen(self, meth, m, n):
        linprog(c = self.c, A_ub=self.A, b_ub=self.b, 
                      method=self.meth)
