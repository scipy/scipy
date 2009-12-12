from numpy.testing import *
import numpy as np
from numpy import array, sqrt
from scipy.special.orthogonal import *
from scipy.special import gamma

from testutils import *

def test_eval_chebyt():
    n = np.arange(0, 10000, 7)
    x = 2*np.random.rand() - 1
    v1 = np.cos(n*np.arccos(x))
    v2 = eval_chebyt(n, x)
    assert np.allclose(v1, v2, rtol=1e-15)

class TestPolys(object):
    """
    Check that the eval_* functions agree with the constructed polynomials

    """

    def check_poly(self, func, cls, param_ranges=[], x_range=[], nn=10,
                   nparam=10, nx=10, rtol=1e-8):
        np.random.seed(1234)

        dataset = []
        for n in np.arange(nn):
            params = [a + (b-a)*np.random.rand(nparam) for a,b in param_ranges]
            params = np.asarray(params).T
            if not param_ranges:
                params = [0]
            for p in params:
                if param_ranges:
                    p = (n,) + tuple(p)
                else:
                    p = (n,)
                x = x_range[0] + (x_range[1] - x_range[0])*np.random.rand(nx)
                poly = np.poly1d(cls(*p))
                z = np.c_[np.tile(p, (nx,1)), x, poly(x)]
                dataset.append(z)

        dataset = np.concatenate(dataset, axis=0)

        def polyfunc(*p):
            p = (p[0].astype(int),) + p[1:]
            return func(*p)

        ds = FuncData(polyfunc, dataset, range(len(param_ranges)+2), -1,
                      rtol=rtol)
        ds.check()

    def test_jacobi(self):
        self.check_poly(eval_jacobi, jacobi,
                   param_ranges=[(-0.99, 10), (-0.99, 10)], x_range=[-1, 1],
                   rtol=1e-5)

    def test_sh_jacobi(self):
        self.check_poly(eval_sh_jacobi, sh_jacobi,
                   param_ranges=[(1, 10), (0, 1)], x_range=[0, 1],
                   rtol=1e-5)

    def test_gegenbauer(self):
        self.check_poly(eval_gegenbauer, gegenbauer,
                   param_ranges=[(-0.499, 10)], x_range=[-1, 1],
                   rtol=1e-7)

    def test_chebyt(self):
        self.check_poly(eval_chebyt, chebyt,
                   param_ranges=[], x_range=[-1, 1])

    def test_chebyu(self):
        self.check_poly(eval_chebyu, chebyu,
                   param_ranges=[], x_range=[-1, 1])

    def test_chebys(self):
        self.check_poly(eval_chebys, chebys,
                   param_ranges=[], x_range=[-2, 2])

    def test_chebyc(self):
        self.check_poly(eval_chebyc, chebyc,
                   param_ranges=[], x_range=[-2, 2])

    def test_sh_chebyt(self):
        self.check_poly(eval_sh_chebyt, sh_chebyt,
                   param_ranges=[], x_range=[0, 1])

    def test_sh_chebyu(self):
        self.check_poly(eval_sh_chebyu, sh_chebyu,
                   param_ranges=[], x_range=[0, 1])

    def test_legendre(self):
        self.check_poly(eval_legendre, legendre,
                   param_ranges=[], x_range=[-1, 1])

    def test_sh_legendre(self):
        self.check_poly(eval_sh_legendre, sh_legendre,
                   param_ranges=[], x_range=[0, 1])

    def test_genlaguerre(self):
        self.check_poly(eval_genlaguerre, genlaguerre,
                   param_ranges=[(-0.99, 10)], x_range=[0, 100])

    def test_laguerre(self):
        self.check_poly(eval_laguerre, laguerre,
                   param_ranges=[], x_range=[0, 100])

    def test_hermite(self):
        self.check_poly(eval_hermite, hermite,
                   param_ranges=[], x_range=[-100, 100])

    def test_hermitenorm(self):
        self.check_poly(eval_hermitenorm, hermitenorm,
                        param_ranges=[], x_range=[-100, 100])
