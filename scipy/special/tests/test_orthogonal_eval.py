import numpy as np
from numpy.testing import assert_
import scipy.special.orthogonal as orth

from scipy.special._testutils import FuncData


def test_eval_chebyt():
    n = np.arange(0, 10000, 7)
    x = 2*np.random.rand() - 1
    v1 = np.cos(n*np.arccos(x))
    v2 = orth.eval_chebyt(n, x)
    assert_(np.allclose(v1, v2, rtol=1e-15))


def test_warnings():
    # ticket 1334
    olderr = np.seterr(all='raise')
    try:
        # these should raise no fp warnings
        orth.eval_legendre(1, 0)
        orth.eval_laguerre(1, 1)
        orth.eval_gegenbauer(1, 1, 0)
    finally:
        np.seterr(**olderr)


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
                x[0] = x_range[0] # always include domain start point
                x[1] = x_range[1] # always include domain end point
                poly = np.poly1d(cls(*p))
                z = np.c_[np.tile(p, (nx,1)), x, poly(x)]
                dataset.append(z)

        dataset = np.concatenate(dataset, axis=0)

        def polyfunc(*p):
            p = (p[0].astype(int),) + p[1:]
            return func(*p)

        olderr = np.seterr(all='raise')
        try:
            ds = FuncData(polyfunc, dataset, range(len(param_ranges)+2), -1,
                          rtol=rtol)
            ds.check()
        finally:
            np.seterr(**olderr)

    def test_jacobi(self):
        self.check_poly(orth.eval_jacobi, orth.jacobi,
                   param_ranges=[(-0.99, 10), (-0.99, 10)], x_range=[-1, 1],
                   rtol=1e-5)

    def test_sh_jacobi(self):
        self.check_poly(orth.eval_sh_jacobi, orth.sh_jacobi,
                   param_ranges=[(1, 10), (0, 1)], x_range=[0, 1],
                   rtol=1e-5)

    def test_gegenbauer(self):
        self.check_poly(orth.eval_gegenbauer, orth.gegenbauer,
                   param_ranges=[(-0.499, 10)], x_range=[-1, 1],
                   rtol=1e-7)

    def test_chebyt(self):
        self.check_poly(orth.eval_chebyt, orth.chebyt,
                   param_ranges=[], x_range=[-1, 1])

    def test_chebyu(self):
        self.check_poly(orth.eval_chebyu, orth.chebyu,
                   param_ranges=[], x_range=[-1, 1])

    def test_chebys(self):
        self.check_poly(orth.eval_chebys, orth.chebys,
                   param_ranges=[], x_range=[-2, 2])

    def test_chebyc(self):
        self.check_poly(orth.eval_chebyc, orth.chebyc,
                   param_ranges=[], x_range=[-2, 2])

    def test_sh_chebyt(self):
        olderr = np.seterr(all='ignore')
        try:
            self.check_poly(orth.eval_sh_chebyt, orth.sh_chebyt,
                            param_ranges=[], x_range=[0, 1])
        finally:
            np.seterr(**olderr)

    def test_sh_chebyu(self):
        self.check_poly(orth.eval_sh_chebyu, orth.sh_chebyu,
                   param_ranges=[], x_range=[0, 1])

    def test_legendre(self):
        self.check_poly(orth.eval_legendre, orth.legendre,
                   param_ranges=[], x_range=[-1, 1])

    def test_sh_legendre(self):
        olderr = np.seterr(all='ignore')
        try:
            self.check_poly(orth.eval_sh_legendre, orth.sh_legendre,
                            param_ranges=[], x_range=[0, 1])
        finally:
            np.seterr(**olderr)

    def test_genlaguerre(self):
        self.check_poly(orth.eval_genlaguerre, orth.genlaguerre,
                   param_ranges=[(-0.99, 10)], x_range=[0, 100])

    def test_laguerre(self):
        self.check_poly(orth.eval_laguerre, orth.laguerre,
                   param_ranges=[], x_range=[0, 100])

    def test_hermite(self):
        self.check_poly(orth.eval_hermite, orth.hermite,
                   param_ranges=[], x_range=[-100, 100])

    def test_hermitenorm(self):
        self.check_poly(orth.eval_hermitenorm, orth.hermitenorm,
                        param_ranges=[], x_range=[-100, 100])
