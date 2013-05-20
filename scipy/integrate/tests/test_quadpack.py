from __future__ import division, print_function, absolute_import

from numpy import sqrt, cos, sin, arctan, exp, log, pi, Inf
from numpy.testing import assert_, TestCase, run_module_suite, dec
from scipy.integrate import quad, dblquad, tplquad, nquad
from scipy.lib.six.moves import xrange
import sys
import math

try:
    import ctypes
    _ctypes_missing = False
except ImportError:
    _ctypes_missing = True


def assert_quad(value_and_err, tabledValue, errTol=1.5e-8):
    value, err = value_and_err
    assert_(abs(value-tabledValue) < err, (value, tabledValue, err))
    if errTol is not None:
        assert_(err < errTol, (err, errTol))


class TestCtypesQuad(TestCase):
    @dec.skipif(_ctypes_missing, msg="Ctypes library could not be found")
    def setUp(self):
        if sys.platform == 'win32':
            file = 'msvcrt.dll'
        elif sys.platform == 'darwin':
            file = 'libm.dylib'
        else:
            file = 'libm.so'
        self.lib = ctypes.CDLL(file)
        restype = ctypes.c_double
        argtypes = (ctypes.c_double,)
        for name in ['sin', 'cos', 'tan']:
            func = getattr(self.lib, name)
            func.restype = restype
            func.argtypes = argtypes

    @dec.skipif(_ctypes_missing, msg="Ctypes library could not be found")
    def test_typical(self):
        assert_quad(quad(self.lib.sin,0,5),quad(math.sin,0,5)[0])
        assert_quad(quad(self.lib.cos,0,5),quad(math.cos,0,5)[0])
        assert_quad(quad(self.lib.tan,0,1),quad(math.tan,0,1)[0])

    #@dec.skipif(_ctypes_missing, msg="Ctypes library could not be found")
    # This doesn't seem to always work.  Need a better way to figure out
    # whether the fast path is called.
    @dec.knownfailureif(True, msg="Unreliable test, see ticket 1684.")
    def test_improvement(self):
        import time
        start = time.time()
        for i in xrange(100):
            quad(self.lib.sin, 0, 100)
        fast = time.time() - start
        start = time.time()
        for i in xrange(100):
            quad(math.sin, 0, 100)
        slow = time.time() - start
        assert_(fast < 0.5*slow, (fast, slow))


class TestQuad(TestCase):
    def test_typical(self):
        # 1) Typical function with two extra arguments:
        def myfunc(x,n,z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        assert_quad(quad(myfunc,0,pi,(2,1.8)), 0.30614353532540296487)

    def test_indefinite(self):
        # 2) Infinite integration limits --- Euler's constant
        def myfunc(x):           # Euler's constant integrand
            return -exp(-x)*log(x)
        assert_quad(quad(myfunc,0,Inf), 0.577215664901532860606512)

    def test_singular(self):
        # 3) Singular points in region of integration.
        def myfunc(x):
            if x > 0 and x < 2.5:
                return sin(x)
            elif x >= 2.5 and x <= 5.0:
                return exp(-x)
            else:
                return 0.0

        assert_quad(quad(myfunc,0,10,points=[2.5,5.0]),
                    1 - cos(2.5) + exp(-2.5) - exp(-5.0))

    def test_sine_weighted_finite(self):
        # 4) Sine weighted integral (finite limits)
        def myfunc(x,a):
            return exp(a*(x-1))

        ome = 2.0**3.4
        assert_quad(quad(myfunc,0,1,args=20,weight='sin',wvar=ome),
                    (20*sin(ome)-ome*cos(ome)+ome*exp(-20))/(20**2 + ome**2))

    def test_sine_weighted_infinite(self):
        # 5) Sine weighted integral (infinite limits)
        def myfunc(x,a):
            return exp(-x*a)

        a = 4.0
        ome = 3.0
        assert_quad(quad(myfunc,0,Inf,args=a,weight='sin',wvar=ome),
                    ome/(a**2 + ome**2))

    def test_cosine_weighted_infinite(self):
        # 6) Cosine weighted integral (negative infinite limits)
        def myfunc(x,a):
            return exp(x*a)

        a = 2.5
        ome = 2.3
        assert_quad(quad(myfunc,-Inf,0,args=a,weight='cos',wvar=ome),
                    a/(a**2 + ome**2))

    def test_algebraic_log_weight(self):
        # 6) Algebraic-logarithmic weight.
        def myfunc(x,a):
            return 1/(1+x+2**(-a))

        a = 1.5
        assert_quad(quad(myfunc,-1,1,args=a,weight='alg',wvar=(-0.5,-0.5)),
                    pi/sqrt((1+2**(-a))**2 - 1))

    def test_cauchypv_weight(self):
        # 7) Cauchy prinicpal value weighting w(x) = 1/(x-c)
        def myfunc(x,a):
            return 2.0**(-a)/((x-1)**2+4.0**(-a))

        a = 0.4
        tabledValue = (2.0**(-0.4)*log(1.5)-2.0**(-1.4)*log((4.0**(-a)+16)/(4.0**(-a)+1))
                     - arctan(2.0**(a+2)) - arctan(2.0**a))/(4.0**(-a) + 1)
        assert_quad(quad(myfunc,0,5,args=0.4,weight='cauchy',wvar=2.0),
                    tabledValue, errTol=1.9e-8)

    def test_double_integral(self):
        # 8) Double Integral test
        def simpfunc(y,x):       # Note order of arguments.
            return x+y

        a, b = 1.0, 2.0
        assert_quad(dblquad(simpfunc,a,b,lambda x: x, lambda x: 2*x),
                    5/6.0 * (b**3.0-a**3.0))

    def test_triple_integral(self):
        # 9) Triple Integral test
        def simpfunc(z,y,x):      # Note order of arguments.
            return x+y+z

        a, b = 1.0, 2.0
        assert_quad(tplquad(simpfunc,a,b,
                            lambda x: x, lambda x: 2*x,
                            lambda x,y: x-y, lambda x,y: x+y),
                    8/3.0 * (b**4.0 - a**4.0))

    def test_n_dimensional_integral(self):
        # 10) Multidimensional Integral test
        def func1(x0,x1,x2,x3):
            return x0**2+x1*x2-x3**3+sin(x0)+(
                1 if (x0-.2*x3-.5-.25*x1>0) else 0)
        def opts_basic(*args):
            return {'points':[.2*args[2]+.5+.25*args[0]]}
        assert_quad(nquad(func1,[[0,1],[-1,1],[.13,.8],[-.15,1]],
                          opts=[opts_basic,{},{},{}]),
                    1.5267454070738635)

        scale = .1
        def func2(x0,x1,x2,x3,t0,t1):
            return x0*x1*x3**2+sin(x2)+1+(1 if x0+t1*x1-t0>0 else 0)
        def lim0(x1,x2,x3,t0,t1):
            return [scale*(x1**2+x2+cos(x3)*t0*t1+1)-1, 
                    scale*(x1**2+x2+cos(x3)*t0*t1+1)+1]
        def lim1(x2,x3,t0,t1):
            return [scale*(t0*x2+t1*x3)-1,
                    scale*(t0*x2+t1*x3)+1]
        def lim2(x3,t0,t1):
            return [scale*(x3+t0**2*t1**3)-1, 
                    scale*(x3+t0**2*t1**3)+1]
        def lim3(t0,t1):
            return [scale*(t0+t1)-1,
                    scale*(t0+t1)+1]
        def opts0(x1,x2,x3,t0,t1):
            return {'points':[t0-t1*x1]}
        def opts1(x2,x3,t0,t1):
            return {}
        def opts2(x3,t0,t1):
            return {}
        def opts3(t0,t1):
            return {}
        assert_quad(nquad(func2,[lim0,lim1,lim2,lim3],args=(0,0),
                          opts=[opts0,opts1,opts2,opts3]),
                    25.066666666666663)


if __name__ == "__main__":
    run_module_suite()
