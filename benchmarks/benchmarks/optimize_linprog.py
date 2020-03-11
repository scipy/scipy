"""
Benchmarks for Linear Programming
"""

# Import testing parameters
from scipy.optimize import linprog, OptimizeWarning
from scipy.linalg import toeplitz
from scipy.optimize.tests.test_linprog import lpgen_2d, magic_square
from numpy.testing import suppress_warnings
from scipy.optimize._remove_redundancy import _remove_redundancy, _remove_redundancy_dense, _remove_redundancy_sparse
from scipy.optimize._linprog_util import _presolve, _clean_inputs, _LPProblem
from scipy.sparse import csc_matrix, csr_matrix, issparse
import numpy as np
import os

from .common import Benchmark

try:
    # the value of SCIPY_XSLOW is used to control whether slow benchmarks run
    slow = int(os.environ.get('SCIPY_XSLOW', 0))
except ValueError:
    pass


methods = [("interior-point", {"sparse": True}),
           ("interior-point", {"sparse": False}),
           ("revised simplex", {})]
rr_methods = [_remove_redundancy, _remove_redundancy_dense,
              _remove_redundancy_sparse]
presolve_methods = ['sparse', 'dense']

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
rr_problems = ['AFIRO', 'BLEND', 'FINNIS', 'RECIPE', 'SCSD6', 'VTP-BASE',
               'BORE3D', 'CYCLE', 'DEGEN2', 'DEGEN3', 'ETAMACRO', 'PILOTNOV',
               'QAP8', 'RECIPE', 'SCORPION', 'SHELL', 'SIERRA', 'WOOD1P']

if not slow:
    problems = ['ADLITTLE', 'AFIRO', 'BLEND', 'BEACONFD', 'GROW7', 'LOTFI',
                'SC105', 'SCTAP1', 'SHARE2B', 'STOCFOR1']
    rr_problems = ['AFIRO', 'BLEND', 'FINNIS', 'RECIPE', 'SCSD6', 'VTP-BASE',
                   'DEGEN2', 'ETAMACRO', 'RECIPE']

presolve_problems = problems


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
    if not slow:
        params[1] = solutions[:2]

    param_names = ['method', '(dimensions, objective)']

    def setup(self, meth, prob):
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


class Netlib_RR(Benchmark):
    params = [
        rr_methods,
        rr_problems
    ]
    param_names = ['method', 'problems']
    # sparse routine returns incorrect matrix on BORE3D and PILOTNOV
    # SVD fails (doesn't converge) on QAP8
    known_fails = {('_remove_redundancy', 'QAP8'),
                   ('_remove_redundancy_sparse', 'BORE3D'),
                   ('_remove_redundancy_sparse', 'PILOTNOV')}

    def setup(self, meth, prob):
        if (meth.__name__, prob) in self.known_fails:
            raise NotImplementedError("Known issues with these benchmarks.")

        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "linprog_benchmark_files",
                                prob + ".npz")
        data = np.load(datafile, allow_pickle=True)

        c, A_eq, A_ub, b_ub, b_eq = (data["c"], data["A_eq"], data["A_ub"],
                                     data["b_ub"], data["b_eq"])
        bounds = np.squeeze(data["bounds"])
        x0 = np.zeros(c.shape)

        lp = _LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
        lp_cleaned = _clean_inputs(lp)
        res = _presolve(lp_cleaned, rr=False, tol=1e-9)[0]

        self.A_eq, self.b_eq = res.A_eq, res.b_eq
        self.true_rank = np.linalg.matrix_rank(self.A_eq)
        if meth == _remove_redundancy_sparse:
            self.A_eq = csc_matrix(self.A_eq)
        self.rr_A = None

    def time_netlib_rr(self, meth, prob):
        self.rr_A, b, status, message = meth(self.A_eq, self.b_eq)

    def track_netlib_rr(self, meth, prob):
        if self.rr_A is None:
            self.time_netlib_rr(meth, prob)

        if meth == _remove_redundancy_sparse:
            self.rr_A = self.rr_A.todense()

        rr_rank = np.linalg.matrix_rank(self.rr_A)
        rr_rows = self.rr_A.shape[0]

        self.error1 = rr_rank - self.true_rank
        self.error2 = rr_rows - self.true_rank

        if abs(self.error1) > abs(self.error2):
            return float(self.error1)
        else:
            return float(self.error2)


class Netlib_presolve(Benchmark):
    params = [
        presolve_methods,
        presolve_problems
    ]
    param_names = ['method', 'problems']

    def setup(self, meth, prob):

        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "linprog_benchmark_files",
                                prob + ".npz")
        data = np.load(datafile, allow_pickle=True)

        c, A_eq, A_ub, b_ub, b_eq = (data["c"], data["A_eq"], data["A_ub"],
                                     data["b_ub"], data["b_eq"])
        bounds = np.squeeze(data["bounds"])
        x0 = np.zeros(c.shape)

        if meth == "sparse":
            A_eq = csr_matrix(A_eq)
            A_ub = csr_matrix(A_ub)

        lp = _LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
        self.lp_cleaned = _clean_inputs(lp)

    def time_netlib_presolve(self, meth, prob):
        _presolve(self.lp_cleaned, rr=False, tol=1e-9)
