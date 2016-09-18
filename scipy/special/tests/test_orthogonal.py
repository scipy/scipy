from __future__ import division, print_function, absolute_import

import numpy as np
from numpy import array, sqrt
from numpy.testing import (assert_array_almost_equal,
                           assert_almost_equal, assert_allclose, assert_raises,
                           run_module_suite)

import scipy.special as sc
from scipy.special.orthogonal import _h_roots_asy
from scipy.special import gamma
from scipy import integrate
from scipy._lib.six import xrange, with_metaclass
from scipy.special._testutils import suppress_warnings, DecoratorMeta


class TestCheby(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_chebyc(self):
        C0 = sc.chebyc(0)
        C1 = sc.chebyc(1)
        olderr = np.seterr(all='ignore')
        try:
            C2 = sc.chebyc(2)
            C3 = sc.chebyc(3)
            C4 = sc.chebyc(4)
            C5 = sc.chebyc(5)
        finally:
            np.seterr(**olderr)

        assert_array_almost_equal(C0.c,[2],13)
        assert_array_almost_equal(C1.c,[1,0],13)
        assert_array_almost_equal(C2.c,[1,0,-2],13)
        assert_array_almost_equal(C3.c,[1,0,-3,0],13)
        assert_array_almost_equal(C4.c,[1,0,-4,0,2],13)
        assert_array_almost_equal(C5.c,[1,0,-5,0,5,0],13)

    def test_chebys(self):
        S0 = sc.chebys(0)
        S1 = sc.chebys(1)
        S2 = sc.chebys(2)
        S3 = sc.chebys(3)
        S4 = sc.chebys(4)
        S5 = sc.chebys(5)
        assert_array_almost_equal(S0.c,[1],13)
        assert_array_almost_equal(S1.c,[1,0],13)
        assert_array_almost_equal(S2.c,[1,0,-1],13)
        assert_array_almost_equal(S3.c,[1,0,-2,0],13)
        assert_array_almost_equal(S4.c,[1,0,-3,0,1],13)
        assert_array_almost_equal(S5.c,[1,0,-4,0,3,0],13)

    def test_chebyt(self):
        T0 = sc.chebyt(0)
        T1 = sc.chebyt(1)
        T2 = sc.chebyt(2)
        T3 = sc.chebyt(3)
        T4 = sc.chebyt(4)
        T5 = sc.chebyt(5)
        assert_array_almost_equal(T0.c,[1],13)
        assert_array_almost_equal(T1.c,[1,0],13)
        assert_array_almost_equal(T2.c,[2,0,-1],13)
        assert_array_almost_equal(T3.c,[4,0,-3,0],13)
        assert_array_almost_equal(T4.c,[8,0,-8,0,1],13)
        assert_array_almost_equal(T5.c,[16,0,-20,0,5,0],13)

    def test_chebyu(self):
        U0 = sc.chebyu(0)
        U1 = sc.chebyu(1)
        U2 = sc.chebyu(2)
        U3 = sc.chebyu(3)
        U4 = sc.chebyu(4)
        U5 = sc.chebyu(5)
        assert_array_almost_equal(U0.c,[1],13)
        assert_array_almost_equal(U1.c,[2,0],13)
        assert_array_almost_equal(U2.c,[4,0,-1],13)
        assert_array_almost_equal(U3.c,[8,0,-4,0],13)
        assert_array_almost_equal(U4.c,[16,0,-12,0,1],13)
        assert_array_almost_equal(U5.c,[32,0,-32,0,6,0],13)


class TestGegenbauer(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_gegenbauer(self):
        a = 5*np.random.random() - 0.5
        if np.any(a == 0):
            a = -0.2
        Ca0 = sc.gegenbauer(0,a)
        Ca1 = sc.gegenbauer(1,a)
        Ca2 = sc.gegenbauer(2,a)
        Ca3 = sc.gegenbauer(3,a)
        Ca4 = sc.gegenbauer(4,a)
        Ca5 = sc.gegenbauer(5,a)

        assert_array_almost_equal(Ca0.c,array([1]),13)
        assert_array_almost_equal(Ca1.c,array([2*a,0]),13)
        assert_array_almost_equal(Ca2.c,array([2*a*(a+1),0,-a]),13)
        assert_array_almost_equal(Ca3.c,array([4*sc.poch(a,3),0,-6*a*(a+1),
                                               0])/3.0,11)
        assert_array_almost_equal(Ca4.c,array([4*sc.poch(a,4),0,-12*sc.poch(a,3),
                                               0,3*a*(a+1)])/6.0,11)
        assert_array_almost_equal(Ca5.c,array([4*sc.poch(a,5),0,-20*sc.poch(a,4),
                                               0,15*sc.poch(a,3),0])/15.0,11)


class TestHermite(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_hermite(self):
        H0 = sc.hermite(0)
        H1 = sc.hermite(1)
        H2 = sc.hermite(2)
        H3 = sc.hermite(3)
        H4 = sc.hermite(4)
        H5 = sc.hermite(5)
        assert_array_almost_equal(H0.c,[1],13)
        assert_array_almost_equal(H1.c,[2,0],13)
        assert_array_almost_equal(H2.c,[4,0,-2],13)
        assert_array_almost_equal(H3.c,[8,0,-12,0],13)
        assert_array_almost_equal(H4.c,[16,0,-48,0,12],12)
        assert_array_almost_equal(H5.c,[32,0,-160,0,120,0],12)

    def test_hermitenorm(self):
        # He_n(x) = 2**(-n/2) H_n(x/sqrt(2))
        psub = np.poly1d([1.0/sqrt(2),0])
        H0 = sc.hermitenorm(0)
        H1 = sc.hermitenorm(1)
        H2 = sc.hermitenorm(2)
        H3 = sc.hermitenorm(3)
        H4 = sc.hermitenorm(4)
        H5 = sc.hermitenorm(5)
        he0 = sc.hermite(0)(psub)
        he1 = sc.hermite(1)(psub) / sqrt(2)
        he2 = sc.hermite(2)(psub) / 2.0
        he3 = sc.hermite(3)(psub) / (2*sqrt(2))
        he4 = sc.hermite(4)(psub) / 4.0
        he5 = sc.hermite(5)(psub) / (4.0*sqrt(2))

        assert_array_almost_equal(H0.c,he0.c,13)
        assert_array_almost_equal(H1.c,he1.c,13)
        assert_array_almost_equal(H2.c,he2.c,13)
        assert_array_almost_equal(H3.c,he3.c,13)
        assert_array_almost_equal(H4.c,he4.c,13)
        assert_array_almost_equal(H5.c,he5.c,13)


class _test_sh_legendre(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_sh_legendre(self):
        # P*_n(x) = P_n(2x-1)
        psub = np.poly1d([2,-1])
        Ps0 = sc.sh_legendre(0)
        Ps1 = sc.sh_legendre(1)
        Ps2 = sc.sh_legendre(2)
        Ps3 = sc.sh_legendre(3)
        Ps4 = sc.sh_legendre(4)
        Ps5 = sc.sh_legendre(5)
        pse0 = sc.legendre(0)(psub)
        pse1 = sc.legendre(1)(psub)
        pse2 = sc.legendre(2)(psub)
        pse3 = sc.legendre(3)(psub)
        pse4 = sc.legendre(4)(psub)
        pse5 = sc.legendre(5)(psub)
        assert_array_almost_equal(Ps0.c,pse0.c,13)
        assert_array_almost_equal(Ps1.c,pse1.c,13)
        assert_array_almost_equal(Ps2.c,pse2.c,13)
        assert_array_almost_equal(Ps3.c,pse3.c,13)
        assert_array_almost_equal(Ps4.c,pse4.c,12)
        assert_array_almost_equal(Ps5.c,pse5.c,12)


class _test_sh_chebyt(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_sh_chebyt(self):
        # T*_n(x) = T_n(2x-1)
        psub = np.poly1d([2,-1])
        Ts0 = sc.sh_chebyt(0)
        Ts1 = sc.sh_chebyt(1)
        Ts2 = sc.sh_chebyt(2)
        Ts3 = sc.sh_chebyt(3)
        Ts4 = sc.sh_chebyt(4)
        Ts5 = sc.sh_chebyt(5)
        tse0 = sc.chebyt(0)(psub)
        tse1 = sc.chebyt(1)(psub)
        tse2 = sc.chebyt(2)(psub)
        tse3 = sc.chebyt(3)(psub)
        tse4 = sc.chebyt(4)(psub)
        tse5 = sc.chebyt(5)(psub)
        assert_array_almost_equal(Ts0.c,tse0.c,13)
        assert_array_almost_equal(Ts1.c,tse1.c,13)
        assert_array_almost_equal(Ts2.c,tse2.c,13)
        assert_array_almost_equal(Ts3.c,tse3.c,13)
        assert_array_almost_equal(Ts4.c,tse4.c,12)
        assert_array_almost_equal(Ts5.c,tse5.c,12)


class _test_sh_chebyu(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_sh_chebyu(self):
        # U*_n(x) = U_n(2x-1)
        psub = np.poly1d([2,-1])
        Us0 = sc.sh_chebyu(0)
        Us1 = sc.sh_chebyu(1)
        Us2 = sc.sh_chebyu(2)
        Us3 = sc.sh_chebyu(3)
        Us4 = sc.sh_chebyu(4)
        Us5 = sc.sh_chebyu(5)
        use0 = sc.chebyu(0)(psub)
        use1 = sc.chebyu(1)(psub)
        use2 = sc.chebyu(2)(psub)
        use3 = sc.chebyu(3)(psub)
        use4 = sc.chebyu(4)(psub)
        use5 = sc.chebyu(5)(psub)
        assert_array_almost_equal(Us0.c,use0.c,13)
        assert_array_almost_equal(Us1.c,use1.c,13)
        assert_array_almost_equal(Us2.c,use2.c,13)
        assert_array_almost_equal(Us3.c,use3.c,13)
        assert_array_almost_equal(Us4.c,use4.c,12)
        assert_array_almost_equal(Us5.c,use5.c,11)


class _test_sh_jacobi(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_sh_jacobi(self):
        # G^(p,q)_n(x) = n! gamma(n+p)/gamma(2*n+p) * P^(p-q,q-1)_n(2*x-1)
        conv = lambda n,p: gamma(n+1)*gamma(n+p)/gamma(2*n+p)
        psub = np.poly1d([2,-1])
        q = 4 * np.random.random()
        p = q-1 + 2*np.random.random()
        #print "shifted jacobi p,q = ", p, q
        G0 = sc.sh_jacobi(0,p,q)
        G1 = sc.sh_jacobi(1,p,q)
        G2 = sc.sh_jacobi(2,p,q)
        G3 = sc.sh_jacobi(3,p,q)
        G4 = sc.sh_jacobi(4,p,q)
        G5 = sc.sh_jacobi(5,p,q)
        ge0 = sc.jacobi(0,p-q,q-1)(psub) * conv(0,p)
        ge1 = sc.jacobi(1,p-q,q-1)(psub) * conv(1,p)
        ge2 = sc.jacobi(2,p-q,q-1)(psub) * conv(2,p)
        ge3 = sc.jacobi(3,p-q,q-1)(psub) * conv(3,p)
        ge4 = sc.jacobi(4,p-q,q-1)(psub) * conv(4,p)
        ge5 = sc.jacobi(5,p-q,q-1)(psub) * conv(5,p)

        assert_array_almost_equal(G0.c,ge0.c,13)
        assert_array_almost_equal(G1.c,ge1.c,13)
        assert_array_almost_equal(G2.c,ge2.c,13)
        assert_array_almost_equal(G3.c,ge3.c,13)
        assert_array_almost_equal(G4.c,ge4.c,13)
        assert_array_almost_equal(G5.c,ge5.c,13)


class TestCall(with_metaclass(DecoratorMeta, object)):
    decorators = [(suppress_warnings, None)]

    def test_call(self):
        poly = []
        for n in xrange(5):
            poly.extend([x.strip() for x in
                ("""
                sc.jacobi(%(n)d,0.3,0.9)
                sc.sh_jacobi(%(n)d,0.3,0.9)
                sc.genlaguerre(%(n)d,0.3)
                sc.laguerre(%(n)d)
                sc.hermite(%(n)d)
                sc.hermitenorm(%(n)d)
                sc.gegenbauer(%(n)d,0.3)
                sc.chebyt(%(n)d)
                sc.chebyu(%(n)d)
                sc.chebyc(%(n)d)
                sc.chebys(%(n)d)
                sc.sh_chebyt(%(n)d)
                sc.sh_chebyu(%(n)d)
                sc.legendre(%(n)d)
                sc.sh_legendre(%(n)d)
                """ % dict(n=n)).split()
            ])
        olderr = np.seterr(all='ignore')
        try:
            for pstr in poly:
                p = eval(pstr)
                assert_almost_equal(p(0.315), np.poly1d(p)(0.315), err_msg=pstr)
        finally:
            np.seterr(**olderr)

def verify_gauss_quad(root_func, eval_func, weight_func, a, b, N,
                      rtol=1e-15, atol=1e-14):
        # this test is copied from numpy's TestGauss in test_hermite.py
        x, w, mu = root_func(N, True)

        n = np.arange(N)
        v = eval_func(n[:,np.newaxis], x)
        vv = np.dot(v*w, v.T)
        vd = 1 / np.sqrt(vv.diagonal())
        vv = vd[:, np.newaxis] * vv * vd
        assert_allclose(vv, np.eye(N), rtol, atol)

        # check that the integral of 1 is correct
        assert_allclose(w.sum(), mu, rtol, atol)

        # compare the results of integrating a function with quad.
        f = lambda x: x**3 - 3*x**2 + x - 2
        resI = integrate.quad(lambda x: f(x)*weight_func(x), a, b)
        resG = np.vdot(f(x), w)
        rtol = 1e-6 if 1e-6 < resI[1] else resI[1] * 10
        assert_allclose(resI[0], resG, rtol=rtol)

def test_j_roots():
    rf = lambda a, b: lambda n, mu: sc.j_roots(n, a, b, mu)
    ef = lambda a, b: lambda n, x: sc.eval_jacobi(n, a, b, x)
    wf = lambda a, b: lambda x: (1 - x)**a * (1 + x)**b

    vgq = verify_gauss_quad
    vgq(rf(-0.5, -0.75), ef(-0.5, -0.75), wf(-0.5, -0.75), -1., 1., 5)
    vgq(rf(-0.5, -0.75), ef(-0.5, -0.75), wf(-0.5, -0.75), -1., 1.,
        25, atol=1e-12)
    vgq(rf(-0.5, -0.75), ef(-0.5, -0.75), wf(-0.5, -0.75), -1., 1.,
        100, atol=1e-11)

    vgq(rf(0.5, -0.5), ef(0.5, -0.5), wf(0.5, -0.5), -1., 1., 5)
    vgq(rf(0.5, -0.5), ef(0.5, -0.5), wf(0.5, -0.5), -1., 1., 25, atol=1.5e-13)
    vgq(rf(0.5, -0.5), ef(0.5, -0.5), wf(0.5, -0.5), -1., 1., 100, atol=1e-12)

    vgq(rf(1, 0.5), ef(1, 0.5), wf(1, 0.5), -1., 1., 5, atol=2e-13)
    vgq(rf(1, 0.5), ef(1, 0.5), wf(1, 0.5), -1., 1., 25, atol=2e-13)
    vgq(rf(1, 0.5), ef(1, 0.5), wf(1, 0.5), -1., 1., 100, atol=1e-12)

    vgq(rf(0.9, 2), ef(0.9, 2), wf(0.9, 2), -1., 1., 5)
    vgq(rf(0.9, 2), ef(0.9, 2), wf(0.9, 2), -1., 1., 25, atol=1e-13)
    vgq(rf(0.9, 2), ef(0.9, 2), wf(0.9, 2), -1., 1., 100, atol=2e-13)

    vgq(rf(18.24, 27.3), ef(18.24, 27.3), wf(18.24, 27.3), -1., 1., 5)
    vgq(rf(18.24, 27.3), ef(18.24, 27.3), wf(18.24, 27.3), -1., 1., 25)
    vgq(rf(18.24, 27.3), ef(18.24, 27.3), wf(18.24, 27.3), -1., 1.,
        100, atol=1e-13)

    vgq(rf(47.1, -0.2), ef(47.1, -0.2), wf(47.1, -0.2), -1., 1., 5, atol=1e-13)
    vgq(rf(47.1, -0.2), ef(47.1, -0.2), wf(47.1, -0.2), -1., 1., 25, atol=2e-13)
    vgq(rf(47.1, -0.2), ef(47.1, -0.2), wf(47.1, -0.2), -1., 1.,
        100, atol=1e-11)

    vgq(rf(2.25, 68.9), ef(2.25, 68.9), wf(2.25, 68.9), -1., 1., 5)
    vgq(rf(2.25, 68.9), ef(2.25, 68.9), wf(2.25, 68.9), -1., 1., 25, atol=1e-13)
    vgq(rf(2.25, 68.9), ef(2.25, 68.9), wf(2.25, 68.9), -1., 1.,
        100, atol=1e-13)

    # when alpha == beta == 0, P_n^{a,b}(x) == P_n(x)
    xj, wj = sc.j_roots(6, 0.0, 0.0)
    xl, wl = sc.p_roots(6)
    assert_allclose(xj, xl, 1e-14, 1e-14)
    assert_allclose(wj, wl, 1e-14, 1e-14)

    # when alpha == beta != 0, P_n^{a,b}(x) == C_n^{alpha+0.5}(x)
    xj, wj = sc.j_roots(6, 4.0, 4.0)
    xc, wc = sc.cg_roots(6, 4.5)
    assert_allclose(xj, xc, 1e-14, 1e-14)
    assert_allclose(wj, wc, 1e-14, 1e-14)

    x, w = sc.j_roots(5, 2, 3, False)
    y, v, m = sc.j_roots(5, 2, 3, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(wf(2,3), -1, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.j_roots, 0, 1, 1)
    assert_raises(ValueError, sc.j_roots, 3.3, 1, 1)
    assert_raises(ValueError, sc.j_roots, 3, -2, 1)
    assert_raises(ValueError, sc.j_roots, 3, 1, -2)
    assert_raises(ValueError, sc.j_roots, 3, -2, -2)

def test_js_roots():
    rf = lambda a, b: lambda n, mu: sc.js_roots(n, a, b, mu)
    ef = lambda a, b: lambda n, x: sc.eval_sh_jacobi(n, a, b, x)
    wf = lambda a, b: lambda x: (1. - x)**(a - b) * (x)**(b - 1.)

    vgq = verify_gauss_quad
    vgq(rf(-0.5, 0.25), ef(-0.5, 0.25), wf(-0.5, 0.25), 0., 1., 5)
    vgq(rf(-0.5, 0.25), ef(-0.5, 0.25), wf(-0.5, 0.25), 0., 1.,
        25, atol=1e-12)
    vgq(rf(-0.5, 0.25), ef(-0.5, 0.25), wf(-0.5, 0.25), 0., 1.,
        100, atol=1e-11)

    vgq(rf(0.5, 0.5), ef(0.5, 0.5), wf(0.5, 0.5), 0., 1., 5)
    vgq(rf(0.5, 0.5), ef(0.5, 0.5), wf(0.5, 0.5), 0., 1., 25, atol=1e-13)
    vgq(rf(0.5, 0.5), ef(0.5, 0.5), wf(0.5, 0.5), 0., 1., 100, atol=1e-12)

    vgq(rf(1, 0.5), ef(1, 0.5), wf(1, 0.5), 0., 1., 5)
    vgq(rf(1, 0.5), ef(1, 0.5), wf(1, 0.5), 0., 1., 25, atol=1.5e-13)
    vgq(rf(1, 0.5), ef(1, 0.5), wf(1, 0.5), 0., 1., 100, atol=1e-12)

    vgq(rf(2, 0.9), ef(2, 0.9), wf(2, 0.9), 0., 1., 5)
    vgq(rf(2, 0.9), ef(2, 0.9), wf(2, 0.9), 0., 1., 25, atol=1e-13)
    vgq(rf(2, 0.9), ef(2, 0.9), wf(2, 0.9), 0., 1., 100, atol=1e-12)

    vgq(rf(27.3, 18.24), ef(27.3, 18.24), wf(27.3, 18.24), 0., 1., 5)
    vgq(rf(27.3, 18.24), ef(27.3, 18.24), wf(27.3, 18.24), 0., 1., 25)
    vgq(rf(27.3, 18.24), ef(27.3, 18.24), wf(27.3, 18.24), 0., 1.,
        100, atol=1e-13)

    vgq(rf(47.1, 0.2), ef(47.1, 0.2), wf(47.1, 0.2), 0., 1., 5, atol=1e-12)
    vgq(rf(47.1, 0.2), ef(47.1, 0.2), wf(47.1, 0.2), 0., 1., 25, atol=1e-11)
    vgq(rf(47.1, 0.2), ef(47.1, 0.2), wf(47.1, 0.2), 0., 1., 100, atol=1e-10)

    vgq(rf(68.9, 2.25), ef(68.9, 2.25), wf(68.9, 2.25), 0., 1., 5, atol=3.5e-14)
    vgq(rf(68.9, 2.25), ef(68.9, 2.25), wf(68.9, 2.25), 0., 1., 25, atol=2e-13)
    vgq(rf(68.9, 2.25), ef(68.9, 2.25), wf(68.9, 2.25), 0., 1.,
        100, atol=1e-12)

    x, w = sc.js_roots(5, 3, 2, False)
    y, v, m = sc.js_roots(5, 3, 2, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(wf(3,2), 0, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.js_roots, 0, 1, 1)
    assert_raises(ValueError, sc.js_roots, 3.3, 1, 1)
    assert_raises(ValueError, sc.js_roots, 3, 1, 2)    # p - q <= -1
    assert_raises(ValueError, sc.js_roots, 3, 2, -1)   # q <= 0
    assert_raises(ValueError, sc.js_roots, 3, -2, -1)  # both

def test_h_roots():
    rootf = sc.h_roots
    evalf = sc.eval_hermite
    weightf = sc.weight_hermite

    verify_gauss_quad(rootf, evalf, weightf, -np.inf, np.inf, 5)
    verify_gauss_quad(rootf, evalf, weightf, -np.inf, np.inf, 25, atol=1e-13)
    verify_gauss_quad(rootf, evalf, weightf, -np.inf, np.inf, 100, atol=1e-12)

    # Golub-Welsch branch
    x, w = sc.h_roots(5, False)
    y, v, m = sc.h_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, -np.inf, np.inf)
    assert_allclose(m, muI, rtol=muI_err)

    # Asymptotic branch (switch over at n >= 150)
    x, w = sc.h_roots(200, False)
    y, v, m = sc.h_roots(200, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)
    assert_allclose(sum(v), m, 1e-14, 1e-14)

    assert_raises(ValueError, sc.h_roots, 0)
    assert_raises(ValueError, sc.h_roots, 3.3)

def test_h_roots_asy():
    # Recursion for Hermite functions
    def hermite_recursion(n, nodes):
        H = np.zeros((n, nodes.size))
        H[0,:] = np.pi**(-0.25) * np.exp(-0.5*nodes**2)
        if n > 1:
            H[1,:] = sqrt(2.0) * nodes * H[0,:]
            for k in xrange(2, n):
                H[k,:] = sqrt(2.0/k) * nodes * H[k-1,:] - sqrt((k-1.0)/k) * H[k-2,:]
        return H

    # This tests only the nodes
    def test(N, rtol=1e-15, atol=1e-14):
        x, w = _h_roots_asy(N)
        H = hermite_recursion(N+1, x)
        assert_allclose(H[-1,:], np.zeros(N), rtol, atol)
        assert_allclose(sum(w), sqrt(np.pi), rtol, atol)

    test(150, atol=1e-12)
    test(151, atol=1e-12)
    test(300, atol=1e-12)
    test(301, atol=1e-12)
    test(500, atol=1e-12)
    test(501, atol=1e-12)
    test(999, atol=1e-12)
    test(1000, atol=1e-12)
    test(2000, atol=1e-12)
    test(5000, atol=1e-12)

def test_he_roots():
    rootf = sc.he_roots
    evalf = sc.eval_hermitenorm
    weightf = sc.weight_hermitenorm

    verify_gauss_quad(rootf, evalf, weightf, -np.inf, np.inf, 5)
    verify_gauss_quad(rootf, evalf, weightf, -np.inf, np.inf, 25, atol=1e-13)
    verify_gauss_quad(rootf, evalf, weightf, -np.inf, np.inf, 100, atol=1e-12)

    x, w = sc.he_roots(5, False)
    y, v, m = sc.he_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, -np.inf, np.inf)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.he_roots, 0)
    assert_raises(ValueError, sc.he_roots, 3.3)

def test_cg_roots():
    rootf = lambda a: lambda n, mu: sc.cg_roots(n, a, mu)
    evalf = lambda a: lambda n, x: sc.eval_gegenbauer(n, a, x)
    weightf = lambda a: lambda x: sc.weight_gegenbauer(a, x)

    vgq = verify_gauss_quad
    vgq(rootf(-0.25), evalf(-0.25), weightf(-0.25), -1., 1., 5)
    vgq(rootf(-0.25), evalf(-0.25), weightf(-0.25), -1., 1., 25, atol=1e-12)
    vgq(rootf(-0.25), evalf(-0.25), weightf(-0.25), -1., 1., 100, atol=1e-11)

    vgq(rootf(0.1), evalf(0.1), weightf(0.1), -1., 1., 5)
    vgq(rootf(0.1), evalf(0.1), weightf(0.1), -1., 1., 25, atol=1e-13)
    vgq(rootf(0.1), evalf(0.1), weightf(0.1), -1., 1., 100, atol=1e-12)

    vgq(rootf(1), evalf(1), weightf(1), -1., 1., 5)
    vgq(rootf(1), evalf(1), weightf(1), -1., 1., 25, atol=1e-13)
    vgq(rootf(1), evalf(1), weightf(1), -1., 1., 100, atol=1e-12)

    vgq(rootf(10), evalf(10), weightf(10), -1., 1., 5)
    vgq(rootf(10), evalf(10), weightf(10), -1., 1., 25, atol=1e-13)
    vgq(rootf(10), evalf(10), weightf(10), -1., 1., 100, atol=1e-12)

    vgq(rootf(50), evalf(50), weightf(50), -1., 1., 5, atol=1e-13)
    vgq(rootf(50), evalf(50), weightf(50), -1., 1., 25, atol=1e-12)
    vgq(rootf(50), evalf(50), weightf(50), -1., 1., 100, atol=1e-11)

    # this is a special case that the old code supported.
    # when alpha = 0, the gegenbauer polynomial is uniformly 0. but it goes
    # to a scaled down copy of T_n(x) there.
    vgq(rootf(0), sc.eval_chebyt, weightf(0), -1., 1., 5)
    vgq(rootf(0), sc.eval_chebyt, weightf(0), -1., 1., 25)
    vgq(rootf(0), sc.eval_chebyt, weightf(0), -1., 1., 100)

    x, w = sc.cg_roots(5, 2, False)
    y, v, m = sc.cg_roots(5, 2, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf(2), -1, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.cg_roots, 0, 2)
    assert_raises(ValueError, sc.cg_roots, 3.3, 2)
    assert_raises(ValueError, sc.cg_roots, 3, -.75)

def test_t_roots():
    weightf = sc.weight_chebyt
    verify_gauss_quad(sc.t_roots, sc.eval_chebyt, weightf, -1., 1., 5)
    verify_gauss_quad(sc.t_roots, sc.eval_chebyt, weightf, -1., 1., 25)
    verify_gauss_quad(sc.t_roots, sc.eval_chebyt, weightf, -1., 1., 100)

    x, w = sc.t_roots(5, False)
    y, v, m = sc.t_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, -1, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.t_roots, 0)
    assert_raises(ValueError, sc.t_roots, 3.3)

def test_u_roots():
    weightf = sc.weight_chebyu
    verify_gauss_quad(sc.u_roots, sc.eval_chebyu, weightf, -1., 1., 5)
    verify_gauss_quad(sc.u_roots, sc.eval_chebyu, weightf, -1., 1., 25)
    verify_gauss_quad(sc.u_roots, sc.eval_chebyu, weightf, -1., 1., 100)

    x, w = sc.u_roots(5, False)
    y, v, m = sc.u_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, -1, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.u_roots, 0)
    assert_raises(ValueError, sc.u_roots, 3.3)

def test_c_roots():
    weightf = sc.weight_chebyc
    verify_gauss_quad(sc.c_roots, sc.eval_chebyc, weightf, -2., 2., 5)
    verify_gauss_quad(sc.c_roots, sc.eval_chebyc, weightf, -2., 2., 25)
    verify_gauss_quad(sc.c_roots, sc.eval_chebyc, weightf, -2., 2., 100)

    x, w = sc.c_roots(5, False)
    y, v, m = sc.c_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, -2, 2)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.c_roots, 0)
    assert_raises(ValueError, sc.c_roots, 3.3)

def test_s_roots():
    weightf = sc.weight_chebys
    verify_gauss_quad(sc.s_roots, sc.eval_chebys, weightf, -2., 2., 5)
    verify_gauss_quad(sc.s_roots, sc.eval_chebys, weightf, -2., 2., 25)
    verify_gauss_quad(sc.s_roots, sc.eval_chebys, weightf, -2., 2., 100)

    x, w = sc.s_roots(5, False)
    y, v, m = sc.s_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, -2, 2)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.s_roots, 0)
    assert_raises(ValueError, sc.s_roots, 3.3)

def test_ts_roots():
    weightf = sc.weight_sh_chebyt
    verify_gauss_quad(sc.ts_roots, sc.eval_sh_chebyt, weightf, 0., 1., 5)
    verify_gauss_quad(sc.ts_roots, sc.eval_sh_chebyt, weightf, 0., 1., 25)
    verify_gauss_quad(sc.ts_roots, sc.eval_sh_chebyt, weightf, 0., 1.,
                      100, atol=1e-13)

    x, w = sc.ts_roots(5, False)
    y, v, m = sc.ts_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, 0, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.ts_roots, 0)
    assert_raises(ValueError, sc.ts_roots, 3.3)

def test_us_roots():
    weightf = sc.weight_sh_chebyu
    verify_gauss_quad(sc.us_roots, sc.eval_sh_chebyu, weightf, 0., 1., 5)
    verify_gauss_quad(sc.us_roots, sc.eval_sh_chebyu, weightf, 0., 1., 25)
    verify_gauss_quad(sc.us_roots, sc.eval_sh_chebyu, weightf, 0., 1.,
                      100, atol=1e-13)

    x, w = sc.us_roots(5, False)
    y, v, m = sc.us_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, 0, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.us_roots, 0)
    assert_raises(ValueError, sc.us_roots, 3.3)

def test_p_roots():
    weightf = sc.weight_legendre
    verify_gauss_quad(sc.p_roots, sc.eval_legendre, weightf, -1., 1., 5)
    verify_gauss_quad(sc.p_roots, sc.eval_legendre, weightf, -1., 1.,
                      25, atol=1e-13)
    verify_gauss_quad(sc.p_roots, sc.eval_legendre, weightf, -1., 1.,
                      100, atol=1e-12)

    x, w = sc.p_roots(5, False)
    y, v, m = sc.p_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, -1, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.p_roots, 0)
    assert_raises(ValueError, sc.p_roots, 3.3)

def test_ps_roots():
    weightf = sc.weight_sh_legendre
    verify_gauss_quad(sc.ps_roots, sc.eval_sh_legendre, weightf, 0., 1., 5)
    verify_gauss_quad(sc.ps_roots, sc.eval_sh_legendre, weightf, 0., 1.,
                      25, atol=1e-13)
    verify_gauss_quad(sc.ps_roots, sc.eval_sh_legendre, weightf, 0., 1.,
                      100, atol=1e-12)

    x, w = sc.ps_roots(5, False)
    y, v, m = sc.ps_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, 0, 1)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.ps_roots, 0)
    assert_raises(ValueError, sc.ps_roots, 3.3)

def test_l_roots():
    weightf = sc.weight_laguerre
    verify_gauss_quad(sc.l_roots, sc.eval_laguerre, weightf, 0., np.inf, 5)
    verify_gauss_quad(sc.l_roots, sc.eval_laguerre, weightf, 0., np.inf,
                      25, atol=1e-13)
    verify_gauss_quad(sc.l_roots, sc.eval_laguerre, weightf, 0., np.inf,
                      100, atol=1e-12)

    x, w = sc.l_roots(5, False)
    y, v, m = sc.l_roots(5, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf, 0, np.inf)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.l_roots, 0)
    assert_raises(ValueError, sc.l_roots, 3.3)

def test_la_roots():
    rootf = lambda a: lambda n, mu: sc.la_roots(n, a, mu)
    evalf = lambda a: lambda n, x: sc.eval_genlaguerre(n, a, x)
    weightf = lambda a: lambda x: sc.weight_genlaguerre(a, x)

    vgq = verify_gauss_quad
    vgq(rootf(-0.5), evalf(-0.5), weightf(-0.5), 0., np.inf, 5)
    vgq(rootf(-0.5), evalf(-0.5), weightf(-0.5), 0., np.inf, 25, atol=1e-13)
    vgq(rootf(-0.5), evalf(-0.5), weightf(-0.5), 0., np.inf, 100, atol=1e-12)

    vgq(rootf(0.1), evalf(0.1), weightf(0.1), 0., np.inf, 5)
    vgq(rootf(0.1), evalf(0.1), weightf(0.1), 0., np.inf, 25, atol=1e-13)
    vgq(rootf(0.1), evalf(0.1), weightf(0.1), 0., np.inf, 100, atol=1e-13)

    vgq(rootf(1), evalf(1), weightf(1), 0., np.inf, 5)
    vgq(rootf(1), evalf(1), weightf(1), 0., np.inf, 25, atol=1e-13)
    vgq(rootf(1), evalf(1), weightf(1), 0., np.inf, 100, atol=1e-13)

    vgq(rootf(10), evalf(10), weightf(10), 0., np.inf, 5)
    vgq(rootf(10), evalf(10), weightf(10), 0., np.inf, 25, atol=1e-13)
    vgq(rootf(10), evalf(10), weightf(10), 0., np.inf, 100, atol=1e-12)

    vgq(rootf(50), evalf(50), weightf(50), 0., np.inf, 5)
    vgq(rootf(50), evalf(50), weightf(50), 0., np.inf, 25, atol=1e-13)
    vgq(rootf(50), evalf(50), weightf(50), 0., np.inf, 100, rtol=1e-14, atol=2e-13)

    x, w = sc.la_roots(5, 2, False)
    y, v, m = sc.la_roots(5, 2, True)
    assert_allclose(x, y, 1e-14, 1e-14)
    assert_allclose(w, v, 1e-14, 1e-14)

    muI, muI_err = integrate.quad(weightf(2.), 0., np.inf)
    assert_allclose(m, muI, rtol=muI_err)

    assert_raises(ValueError, sc.la_roots, 0, 2)
    assert_raises(ValueError, sc.la_roots, 3.3, 2)
    assert_raises(ValueError, sc.la_roots, 3, -1.1)


if __name__ == "__main__":
    run_module_suite()
