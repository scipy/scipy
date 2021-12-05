import os

import numpy as np
from numpy.testing import suppress_warnings

from .common import Benchmark, is_xslow, safe_import

with safe_import():
    from scipy.optimize import linprog, OptimizeWarning

with safe_import():
    from scipy.optimize.tests.test_linprog import lpgen_2d, magic_square

with safe_import():
    from scipy.linalg import toeplitz

methods = [("highs-ipm", {}),
           ("highs-ds", {})]

problems = ['25FV47', '80BAU3B', 'ADLITTLE', 'AFIRO', 'AGG', 'AGG2', 'AGG3',
            'BANDM', 'BEACONFD', 'BLEND', 'BNL1', 'BNL2', 'BORE3D', 'BRANDY',
            'CAPRI', 'CYCLE', 'CZPROB', 'D2Q06C', 'D6CUBE', 'DEGEN2', 'DEGEN3',
            'DFL001', 'E226', 'ETAMACRO', 'FFFFF800', 'FINNIS', 'FIT1D',
            'FIT1P', 'FIT2D', 'FIT2P', 'GANGES', 'GFRD-PNC', 'GREENBEA',
            'GREENBEB', 'GROW15', 'GROW22', 'GROW7', 'ISRAEL', 'KB2', 'LOTFI',
            'MAROS', 'MAROS-R7', 'MODSZK1', 'PEROLD', 'PILOT', 'PILOT4',
            'PILOT87', 'PILOT-JA', 'PILOTNOV', 'PILOT-WE', 'QAP8', 'QAP12',
            'QAP15', 'RECIPE', 'SC105', 'SC205', 'SC50A', 'SC50B', 'SCAGR25',
            'SCAGR7', 'SCFXM1', 'SCFXM2', 'SCFXM3', 'SCORPION', 'SCRS8',
            'SCSD1', 'SCSD6', 'SCSD8', 'SCTAP1', 'SCTAP2', 'SCTAP3', 'SHARE1B',
            'SHARE2B', 'SHELL', 'SHIP04L', 'SHIP04S', 'SHIP08L', 'SHIP08S',
            'SHIP12L', 'SHIP12S', 'SIERRA', 'STAIR', 'STANDATA', 'STANDMPS',
            'STOCFOR1', 'STOCFOR2', 'STOCFOR3', 'TRUSS', 'TUFF', 'VTP-BASE',
            'WOOD1P', 'WOODW']
infeasible_problems = ['bgdbg1', 'bgetam', 'bgindy', 'bgprtr', 'box1',
                       'ceria3d', 'chemcom', 'cplex1', 'cplex2', 'ex72a',
                       'ex73a', 'forest6', 'galenet', 'gosh', 'gran',
                       'itest2', 'itest6', 'klein1', 'klein2', 'klein3',
                       'mondou2', 'pang', 'pilot4i', 'qual', 'reactor',
                       'refinery', 'vol1', 'woodinfe']

if not is_xslow():
    enabled_problems = ['ADLITTLE', 'AFIRO', 'BLEND', 'BEACONFD', 'GROW7',
                        'LOTFI', 'SC105', 'SCTAP1', 'SHARE2B', 'STOCFOR1']
    enabled_infeasible_problems = ['bgdbg1', 'bgprtr', 'box1', 'chemcom',
                                   'cplex2', 'ex72a', 'ex73a', 'forest6',
                                   'galenet', 'itest2', 'itest6', 'klein1',
                                   'refinery', 'woodinfe']
else:
    enabled_problems = problems
    enabled_infeasible_problems = infeasible_problems


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

    solutions = [(3, 1.7305505947214375), (4, 1.5485271031586025),
                 (5, 1.807494583582637), (6, 1.747266446858304)]

    params = [methods, solutions]
    param_names = ['method', '(dimensions, objective)']

    def setup(self, meth, prob):
        if not is_xslow():
            if prob[0] > 4:
                raise NotImplementedError("skipped")

        dims, obj = prob
        self.A_eq, self.b_eq, self.c, numbers = magic_square(dims)
        self.fun = None

    def time_magic_square(self, meth, prob):
        method, options = meth
        with suppress_warnings() as sup:
            sup.filter(OptimizeWarning, "A_eq does not appear")
            res = linprog(c=self.c, A_eq=self.A_eq, b_eq=self.b_eq,
                          bounds=(0, 1), method=method, options=options)
            self.fun = res.fun

    def track_magic_square(self, meth, prob):
        dims, obj = prob
        if self.fun is None:
            self.time_magic_square(meth, prob)
        self.abs_error = np.abs(self.fun - obj)
        self.rel_error = np.abs((self.fun - obj)/obj)
        return min(self.abs_error, self.rel_error)


class KleeMinty(Benchmark):

    params = [
        methods,
        [3, 6, 9]
    ]
    param_names = ['method', 'dimensions']

    def setup(self, meth, dims):
        self.c, self.A_ub, self.b_ub, self.xf, self.obj = klee_minty(dims)
        self.fun = None

    def time_klee_minty(self, meth, dims):
        method, options = meth
        res = linprog(c=self.c, A_ub=self.A_ub, b_ub=self.b_ub,
                      method=method, options=options)
        self.fun = res.fun
        self.x = res.x

    def track_klee_minty(self, meth, prob):
        if self.fun is None:
            self.time_klee_minty(meth, prob)
        self.abs_error = np.abs(self.fun - self.obj)
        self.rel_error = np.abs((self.fun - self.obj)/self.obj)
        return min(self.abs_error, self.rel_error)


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
        if prob not in enabled_problems:
            raise NotImplementedError("skipped")

        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "linprog_benchmark_files",
                                prob + ".npz")
        data = np.load(datafile, allow_pickle=True)
        self.c = data["c"]
        self.A_eq = data["A_eq"]
        self.A_ub = data["A_ub"]
        self.b_ub = data["b_ub"]
        self.b_eq = data["b_eq"]
        self.bounds = np.squeeze(data["bounds"])
        self.obj = float(data["obj"].flatten()[0])
        self.fun = None

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

    def track_netlib(self, meth, prob):
        if self.fun is None:
            self.time_netlib(meth, prob)
        self.abs_error = np.abs(self.fun - self.obj)
        self.rel_error = np.abs((self.fun - self.obj)/self.obj)
        return min(self.abs_error, self.rel_error)


class Netlib_infeasible(Benchmark):
    params = [
        methods,
        infeasible_problems
    ]
    param_names = ['method', 'problems']

    def setup(self, meth, prob):
        if prob not in enabled_infeasible_problems:
            raise NotImplementedError("skipped")

        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "linprog_benchmark_files",
                                "infeasible", prob + ".npz")
        data = np.load(datafile, allow_pickle=True)
        self.c = data["c"]
        self.A_eq = data["A_eq"]
        self.A_ub = data["A_ub"]
        self.b_ub = data["b_ub"]
        self.b_eq = data["b_eq"]
        self.bounds = np.squeeze(data["bounds"])
        self.status = None

    def time_netlib_infeasible(self, meth, prob):
        method, options = meth
        res = linprog(c=self.c,
                      A_ub=self.A_ub,
                      b_ub=self.b_ub,
                      A_eq=self.A_eq,
                      b_eq=self.b_eq,
                      bounds=self.bounds,
                      method=method,
                      options=options)
        self.status = res.status

    def track_netlib_infeasible(self, meth, prob):
        if self.status is None:
            self.time_netlib_infeasible(meth, prob)
        return self.status
