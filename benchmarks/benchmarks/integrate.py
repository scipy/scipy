import numpy as np
from .common import Benchmark, safe_import

from scipy.integrate import quad

with safe_import():
    import ctypes
    import scipy.integrate._test_multivariate as clib_test
    from scipy._lib import _ccallback_c

with safe_import() as exc:
    from scipy import LowLevelCallable
    from_cython = LowLevelCallable.from_cython
if exc.error:
    LowLevelCallable = lambda func, data: (func, data)
    from_cython = lambda *a: a

with safe_import() as exc:
    import cffi
if exc.error:
    cffi = None

with safe_import():
    from scipy.integrate import solve_bvp


class SolveBVP(Benchmark):
    TOL = 1e-5

    def fun_flow(self, x, y, p):
        A = p[0]
        return np.vstack((
            y[1], y[2], 100 * (y[1] ** 2 - y[0] * y[2] - A),
            y[4], -100 * y[0] * y[4] - 1, y[6], -70 * y[0] * y[6]
        ))

    def bc_flow(self, ya, yb, p):
        return np.array([
            ya[0], ya[1], yb[0] - 1, yb[1], ya[3], yb[3], ya[5], yb[5] - 1])

    def time_flow(self):
        x = np.linspace(0, 1, 10)
        y = np.ones((7, x.size))
        solve_bvp(self.fun_flow, self.bc_flow, x, y, p=[1], tol=self.TOL)

    def fun_peak(self, x, y):
        eps = 1e-3
        return np.vstack((
            y[1],
            -(4 * x * y[1] + 2 * y[0]) / (eps + x**2)
        ))

    def bc_peak(self, ya, yb):
        eps = 1e-3
        v = (1 + eps) ** -1
        return np.array([ya[0] - v, yb[0] - v])

    def time_peak(self):
        x = np.linspace(-1, 1, 5)
        y = np.zeros((2, x.size))
        solve_bvp(self.fun_peak, self.bc_peak, x, y, tol=self.TOL)

    def fun_gas(self, x, y):
        alpha = 0.8
        return np.vstack((
            y[1],
            -2 * x * y[1] * (1 - alpha * y[0]) ** -0.5
        ))

    def bc_gas(self, ya, yb):
        return np.array([ya[0] - 1, yb[0]])

    def time_gas(self):
        x = np.linspace(0, 3, 5)
        y = np.empty((2, x.size))
        y[0] = 0.5
        y[1] = -0.5
        solve_bvp(self.fun_gas, self.bc_gas, x, y, tol=self.TOL)


class Quad(Benchmark):
    def setup(self):
        from math import sin

        self.f_python = lambda x: sin(x)
        self.f_cython = from_cython(_ccallback_c, "sine")

        try:
            from scipy.integrate.tests.test_quadpack import get_clib_test_routine
            self.f_ctypes = get_clib_test_routine('_multivariate_sin', ctypes.c_double,
                                                  ctypes.c_int, ctypes.c_double)
        except ImportError:
            lib = ctypes.CDLL(clib_test.__file__)
            self.f_ctypes = lib._multivariate_sin
            self.f_ctypes.restype = ctypes.c_double
            self.f_ctypes.argtypes = (ctypes.c_int, ctypes.c_double)

        if cffi is not None:
            voidp = ctypes.cast(self.f_ctypes, ctypes.c_void_p)
            address = voidp.value
            ffi = cffi.FFI()
            self.f_cffi = LowLevelCallable(ffi.cast("double (*)(int, double *)", address))

    def time_quad_python(self):
        quad(self.f_python, 0, np.pi)

    def time_quad_cython(self):
        quad(self.f_cython, 0, np.pi)

    def time_quad_ctypes(self):
        quad(self.f_ctypes, 0, np.pi)

    def time_quad_cffi(self):
        quad(self.f_cffi, 0, np.pi)
