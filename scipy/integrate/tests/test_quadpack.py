from numpy import sqrt, cos, sin, arctan, exp, log, pi, Inf
from numpy.testing import assert_, TestCase, run_module_suite
from scipy.integrate import quad, dblquad, tplquad

def assert_quad((value, err), tabledValue, errTol=1.5e-8):
    assert_(abs(value-tabledValue) < err, (value, tabledValue, err))
    if errTol is not None:
        assert_(err < errTol, (err, errTol))

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
            elif x>= 2.5 and x <= 5.0:
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

if __name__ == "__main__":
    run_module_suite()
