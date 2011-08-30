from numpy.testing import assert_array_almost_equal, assert_almost_equal, \
        rand, TestCase
import numpy as np
from numpy import array, sqrt
import scipy.special.orthogonal as orth
from scipy.special import gamma


class TestCheby(TestCase):
    def test_chebyc(self):
        C0 = orth.chebyc(0)
        C1 = orth.chebyc(1)
        olderr = np.seterr(all='ignore')
        try:
            C2 = orth.chebyc(2)
            C3 = orth.chebyc(3)
            C4 = orth.chebyc(4)
            C5 = orth.chebyc(5)
        finally:
            np.seterr(**olderr)

        assert_array_almost_equal(C0.c,[2],13)
        assert_array_almost_equal(C1.c,[1,0],13)
        assert_array_almost_equal(C2.c,[1,0,-2],13)
        assert_array_almost_equal(C3.c,[1,0,-3,0],13)
        assert_array_almost_equal(C4.c,[1,0,-4,0,2],13)
        assert_array_almost_equal(C5.c,[1,0,-5,0,5,0],13)

    def test_chebys(self):
        S0 = orth.chebys(0)
        S1 = orth.chebys(1)
        S2 = orth.chebys(2)
        S3 = orth.chebys(3)
        S4 = orth.chebys(4)
        S5 = orth.chebys(5)
        assert_array_almost_equal(S0.c,[1],13)
        assert_array_almost_equal(S1.c,[1,0],13)
        assert_array_almost_equal(S2.c,[1,0,-1],13)
        assert_array_almost_equal(S3.c,[1,0,-2,0],13)
        assert_array_almost_equal(S4.c,[1,0,-3,0,1],13)
        assert_array_almost_equal(S5.c,[1,0,-4,0,3,0],13)

    def test_chebyt(self):
        T0 = orth.chebyt(0)
        T1 = orth.chebyt(1)
        T2 = orth.chebyt(2)
        T3 = orth.chebyt(3)
        T4 = orth.chebyt(4)
        T5 = orth.chebyt(5)
        assert_array_almost_equal(T0.c,[1],13)
        assert_array_almost_equal(T1.c,[1,0],13)
        assert_array_almost_equal(T2.c,[2,0,-1],13)
        assert_array_almost_equal(T3.c,[4,0,-3,0],13)
        assert_array_almost_equal(T4.c,[8,0,-8,0,1],13)
        assert_array_almost_equal(T5.c,[16,0,-20,0,5,0],13)

    def test_chebyu(self):
        U0 = orth.chebyu(0)
        U1 = orth.chebyu(1)
        U2 = orth.chebyu(2)
        U3 = orth.chebyu(3)
        U4 = orth.chebyu(4)
        U5 = orth.chebyu(5)
        assert_array_almost_equal(U0.c,[1],13)
        assert_array_almost_equal(U1.c,[2,0],13)
        assert_array_almost_equal(U2.c,[4,0,-1],13)
        assert_array_almost_equal(U3.c,[8,0,-4,0],13)
        assert_array_almost_equal(U4.c,[16,0,-12,0,1],13)
        assert_array_almost_equal(U5.c,[32,0,-32,0,6,0],13)

class TestGegenbauer(TestCase):

    def test_gegenbauer(self):
        a = 5*rand()-0.5
        if np.any(a==0): a = -0.2
        Ca0 = orth.gegenbauer(0,a)
        Ca1 = orth.gegenbauer(1,a)
        Ca2 = orth.gegenbauer(2,a)
        Ca3 = orth.gegenbauer(3,a)
        Ca4 = orth.gegenbauer(4,a)
        Ca5 = orth.gegenbauer(5,a)

        assert_array_almost_equal(Ca0.c,array([1]),13)
        assert_array_almost_equal(Ca1.c,array([2*a,0]),13)
        assert_array_almost_equal(Ca2.c,array([2*a*(a+1),0,-a]),13)
        assert_array_almost_equal(Ca3.c,array([4*orth.poch(a,3),0,-6*a*(a+1),
                                               0])/3.0,11)
        assert_array_almost_equal(Ca4.c,array([4*orth.poch(a,4),0,-12*orth.poch(a,3),
                                               0,3*a*(a+1)])/6.0,11)
        assert_array_almost_equal(Ca5.c,array([4*orth.poch(a,5),0,-20*orth.poch(a,4),
                                               0,15*orth.poch(a,3),0])/15.0,11)

class TestHermite(TestCase):
    def test_hermite(self):
        H0 = orth.hermite(0)
        H1 = orth.hermite(1)
        H2 = orth.hermite(2)
        H3 = orth.hermite(3)
        H4 = orth.hermite(4)
        H5 = orth.hermite(5)
        assert_array_almost_equal(H0.c,[1],13)
        assert_array_almost_equal(H1.c,[2,0],13)
        assert_array_almost_equal(H2.c,[4,0,-2],13)
        assert_array_almost_equal(H3.c,[8,0,-12,0],13)
        assert_array_almost_equal(H4.c,[16,0,-48,0,12],12)
        assert_array_almost_equal(H5.c,[32,0,-160,0,120,0],12)

    def test_hermitenorm(self):
        # He_n(x) = 2**(-n/2) H_n(x/sqrt(2))
        psub = np.poly1d([1.0/sqrt(2),0])
        H0 = orth.hermitenorm(0)
        H1 = orth.hermitenorm(1)
        H2 = orth.hermitenorm(2)
        H3 = orth.hermitenorm(3)
        H4 = orth.hermitenorm(4)
        H5 = orth.hermitenorm(5)
        he0 = orth.hermite(0)(psub)
        he1 = orth.hermite(1)(psub) / sqrt(2)
        he2 = orth.hermite(2)(psub) / 2.0
        he3 = orth.hermite(3)(psub) / (2*sqrt(2))
        he4 = orth.hermite(4)(psub) / 4.0
        he5 = orth.hermite(5)(psub) / (4.0*sqrt(2))

        assert_array_almost_equal(H0.c,he0.c,13)
        assert_array_almost_equal(H1.c,he1.c,13)
        assert_array_almost_equal(H2.c,he2.c,13)
        assert_array_almost_equal(H3.c,he3.c,13)
        assert_array_almost_equal(H4.c,he4.c,13)
        assert_array_almost_equal(H5.c,he5.c,13)

class _test_sh_legendre(TestCase):

    def test_sh_legendre(self):
        # P*_n(x) = P_n(2x-1)
        psub = np.poly1d([2,-1])
        Ps0 = orth.sh_legendre(0)
        Ps1 = orth.sh_legendre(1)
        Ps2 = orth.sh_legendre(2)
        Ps3 = orth.sh_legendre(3)
        Ps4 = orth.sh_legendre(4)
        Ps5 = orth.sh_legendre(5)
        pse0 = orth.legendre(0)(psub)
        pse1 = orth.legendre(1)(psub)
        pse2 = orth.legendre(2)(psub)
        pse3 = orth.legendre(3)(psub)
        pse4 = orth.legendre(4)(psub)
        pse5 = orth.legendre(5)(psub)
        assert_array_almost_equal(Ps0.c,pse0.c,13)
        assert_array_almost_equal(Ps1.c,pse1.c,13)
        assert_array_almost_equal(Ps2.c,pse2.c,13)
        assert_array_almost_equal(Ps3.c,pse3.c,13)
        assert_array_almost_equal(Ps4.c,pse4.c,12)
        assert_array_almost_equal(Ps5.c,pse5.c,12)

class _test_sh_chebyt(TestCase):

    def test_sh_chebyt(self):
        # T*_n(x) = T_n(2x-1)
        psub = np.poly1d([2,-1])
        Ts0 = orth.sh_chebyt(0)
        Ts1 = orth.sh_chebyt(1)
        Ts2 = orth.sh_chebyt(2)
        Ts3 = orth.sh_chebyt(3)
        Ts4 = orth.sh_chebyt(4)
        Ts5 = orth.sh_chebyt(5)
        tse0 = orth.chebyt(0)(psub)
        tse1 = orth.chebyt(1)(psub)
        tse2 = orth.chebyt(2)(psub)
        tse3 = orth.chebyt(3)(psub)
        tse4 = orth.chebyt(4)(psub)
        tse5 = orth.chebyt(5)(psub)
        assert_array_almost_equal(Ts0.c,tse0.c,13)
        assert_array_almost_equal(Ts1.c,tse1.c,13)
        assert_array_almost_equal(Ts2.c,tse2.c,13)
        assert_array_almost_equal(Ts3.c,tse3.c,13)
        assert_array_almost_equal(Ts4.c,tse4.c,12)
        assert_array_almost_equal(Ts5.c,tse5.c,12)

class _test_sh_chebyu(TestCase):

    def test_sh_chebyu(self):
        # U*_n(x) = U_n(2x-1)
        psub = np.poly1d([2,-1])
        Us0 = orth.sh_chebyu(0)
        Us1 = orth.sh_chebyu(1)
        Us2 = orth.sh_chebyu(2)
        Us3 = orth.sh_chebyu(3)
        Us4 = orth.sh_chebyu(4)
        Us5 = orth.sh_chebyu(5)
        use0 = orth.chebyu(0)(psub)
        use1 = orth.chebyu(1)(psub)
        use2 = orth.chebyu(2)(psub)
        use3 = orth.chebyu(3)(psub)
        use4 = orth.chebyu(4)(psub)
        use5 = orth.chebyu(5)(psub)
        assert_array_almost_equal(Us0.c,use0.c,13)
        assert_array_almost_equal(Us1.c,use1.c,13)
        assert_array_almost_equal(Us2.c,use2.c,13)
        assert_array_almost_equal(Us3.c,use3.c,13)
        assert_array_almost_equal(Us4.c,use4.c,12)
        assert_array_almost_equal(Us5.c,use5.c,11)

class _test_sh_jacobi(TestCase):
    def test_sh_jacobi(self):
        # G^(p,q)_n(x) = n! gamma(n+p)/gamma(2*n+p) * P^(p-q,q-1)_n(2*x-1)
        conv = lambda n,p: gamma(n+1)*gamma(n+p)/gamma(2*n+p)
        psub = np.poly1d([2,-1])
        q = 4*rand()
        p = q-1 + 2*rand()
        #print "shifted jacobi p,q = ", p, q
        G0 = orth.sh_jacobi(0,p,q)
        G1 = orth.sh_jacobi(1,p,q)
        G2 = orth.sh_jacobi(2,p,q)
        G3 = orth.sh_jacobi(3,p,q)
        G4 = orth.sh_jacobi(4,p,q)
        G5 = orth.sh_jacobi(5,p,q)
        ge0 = orth.jacobi(0,p-q,q-1)(psub) * conv(0,p)
        ge1 = orth.jacobi(1,p-q,q-1)(psub) * conv(1,p)
        ge2 = orth.jacobi(2,p-q,q-1)(psub) * conv(2,p)
        ge3 = orth.jacobi(3,p-q,q-1)(psub) * conv(3,p)
        ge4 = orth.jacobi(4,p-q,q-1)(psub) * conv(4,p)
        ge5 = orth.jacobi(5,p-q,q-1)(psub) * conv(5,p)

        assert_array_almost_equal(G0.c,ge0.c,13)
        assert_array_almost_equal(G1.c,ge1.c,13)
        assert_array_almost_equal(G2.c,ge2.c,13)
        assert_array_almost_equal(G3.c,ge3.c,13)
        assert_array_almost_equal(G4.c,ge4.c,13)
        assert_array_almost_equal(G5.c,ge5.c,13)

class TestCall(object):
    def test_call(self):
        poly = []
        for n in xrange(5):
            poly.extend([x.strip() for x in
                ("""
                orth.jacobi(%(n)d,0.3,0.9)
                orth.sh_jacobi(%(n)d,0.3,0.9)
                orth.genlaguerre(%(n)d,0.3)
                orth.laguerre(%(n)d)
                orth.hermite(%(n)d)
                orth.hermitenorm(%(n)d)
                orth.gegenbauer(%(n)d,0.3)
                orth.chebyt(%(n)d)
                orth.chebyu(%(n)d)
                orth.chebyc(%(n)d)
                orth.chebys(%(n)d)
                orth.sh_chebyt(%(n)d)
                orth.sh_chebyu(%(n)d)
                orth.legendre(%(n)d)
                orth.sh_legendre(%(n)d)
                """ % dict(n=n)).split()
            ])
        olderr = np.seterr(all='ignore')
        try:
            for pstr in poly:
                p = eval(pstr)
                assert_almost_equal(p(0.315), np.poly1d(p)(0.315), err_msg=pstr)
        finally:
            np.seterr(**olderr)
