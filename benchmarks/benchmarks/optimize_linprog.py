"""
Benchmarks for Linear Programming
"""
from __future__ import division, print_function, absolute_import

# Import testing parameters
try:
    from scipy.optimize import linprog, OptimizeWarning
    from scipy.linalg import toeplitz
    from scipy.optimize.tests.test_linprog import lpgen_2d, magic_square
    from numpy.testing import suppress_warnings
    import numpy as np
    import os
except ImportError:
    pass

from .common import Benchmark

methods = [("interior-point", {"sparse": True}),
           ("interior-point", {"sparse": False}),
           ("revised simplex", {})]
problems = ['25FV47', '80BAU3B', 'ADLITTLE', 'AFIRO', 'AGG', 'AGG2', 'AGG3',
            'BANDM', 'BEACONFD', 'BLEND', 'BNL1', 'BNL2', 'BORE3D', 'BRANDY',
            'CAPRI', 'CYCLE', 'CZPROB', 'D6CUBE', 'DEGEN2', 'DEGEN3', 'E226',
            'ETAMACRO', 'FFFFF800', 'FINNIS', 'FIT1D', 'FIT1P', 'GANGES',
            'GFRD-PNC', 'GROW15', 'GROW22', 'GROW7', 'ISRAEL', 'KB2', 'LOTFI',
            'MAROS', 'MODSZK1', 'PEROLD', 'PILOT', 'PILOT-WE', 'PILOT4',
            'PILOTNOV', 'QAP8', 'RECIPE', 'SC105', 'SC205', 'SC50A', 'SC50B',
            'SCAGR25', 'SCAGR7', 'SCFXM1', 'SCFXM2', 'SCFXM3', 'SCORPION',
            'SCRS8', 'SCSD1', 'SCSD6', 'SCSD8', 'SCTAP1', 'SCTAP2', 'SCTAP3',
            'SHARE1B', 'SHARE2B', 'SHELL', 'SHIP04L', 'SHIP04S', 'SHIP08L',
            'SHIP08S', 'SHIP12L', 'SHIP12S', 'SIERRA', 'STAIR', 'STANDATA',
            'STANDMPS', 'STOCFOR1', 'STOCFOR2', 'TRUSS', 'TUFF', 'VTP-BASE',
            'WOOD1P', 'WOODW']


def klee_minty(D):
    A_1 = np.array([2**(i + 1) if i > 0 else 1 for i in range(D)])
    A1_ = np.zeros(D)
    A1_[0] = 1
    A_ub = toeplitz(A_1, A1_)
    b_ub = np.array([5**(i + 1) for i in range(D)])
    c = -np.array([2**(D - i - 1) for i in range(D)])
    xf = np.zeros(D)
    xf[-1] = 5**D
    obj = c @ xf
    return c, A_ub, b_ub, xf, obj


class MagicSquare(Benchmark):

    params = [
        methods,
        [(3, 1.7305505947214375), (4, 1.5485271031586025)]
    ]
    param_names = ['method', '(dimensions, objective)']

    def setup(self, meth, prob):
        dims, obj = prob
        self.A_eq, self.b_eq, self.c, numbers = magic_square(dims)

    def time_magic_square(self, meth, prob):
        method, options = meth
        with suppress_warnings() as sup:
            sup.filter(OptimizeWarning, "A_eq does not appear")
            res = linprog(c=self.c, A_eq=self.A_eq, b_eq=self.b_eq,
                          bounds=(0, 1), method=method, options=options)
            self.fun = res.fun

    def teardown(self, meth, prob):
        dims, obj = prob
        np.testing.assert_allclose(obj, self.fun, rtol=1e-6, atol=1e-3)


class KleeMinty(Benchmark):

    params = [
        methods,
        [3, 6, 9]
    ]
    param_names = ['method', 'dimensions']

    def setup(self, meth, dims):
        self.c, self.A_ub, self.b_ub, self.xf, self.obj = klee_minty(dims)

    def time_klee_minty(self, meth, dims):
        method, options = meth
        res = linprog(c=self.c, A_ub=self.A_ub, b_ub=self.b_ub,
                      method=method, options=options)
        self.fun = res.fun
        self.x = res.x

    def teardown(self, meth, prob):
        np.testing.assert_allclose(self.obj, self.fun, rtol=1e-6, atol=1e-3)
        np.testing.assert_allclose(self.xf, self.x, rtol=1e-6, atol=1e-3)


class LpGen(Benchmark):
    params = [
        methods,
        range(20, 100, 20),
        range(20, 100, 20)
    ]
    param_names = ['method', 'm', 'n']

    def setup(self, meth, m, n):
        self.A, self.b, self.c = lpgen_2d(m, n)

    def time_lpgen(self, meth, m, n):
        method, options = meth
        with suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "scipy.linalg.solve\nIll-conditioned")
            linprog(c=self.c, A_ub=self.A, b_ub=self.b,
                    method=method, options=options)


class Netlib(Benchmark):
    params = [
        methods,
        problems
    ]
    param_names = ['method', 'problems']

    def setup(self, meth, prob):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        data = np.load(dir_path + "/linprog_benchmark_files/" + prob + ".npz",
                       allow_pickle=True)
        self.c = data["c"]
        self.A_eq = data["A_eq"]
        self.A_ub = data["A_ub"]
        self.b_ub = data["b_ub"]
        self.b_eq = data["b_eq"]
        self.bounds = np.squeeze(data["bounds"])
        self.obj = float(data["obj"].flatten()[0])

    def time_netlib(self, meth, prob):
        method, options = meth
        res = linprog(c=self.c,
                      A_ub=self.A_ub,
                      b_ub=self.b_ub,
                      A_eq=self.A_eq,
                      b_eq=self.b_eq,
                      bounds=self.bounds,
                      method=method,
                      options=options)
        self.fun = res.fun

    def teardown(self, meth, prob):
        np.testing.assert_allclose(self.obj, self.fun, rtol=1e-4, atol=1e-4)
