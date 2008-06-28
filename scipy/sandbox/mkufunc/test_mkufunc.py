import math
import unittest

from numpy import array, arange, allclose, vectorize

from mkufunc import Cfunc, genufunc, mkufunc

class Internal_Tests(unittest.TestCase):
    
    def test_Cfunc(self):
        def sqr(x):
            return x * x
        cf = Cfunc(sqr, [int, int], 42)
        self.assertEqual(cf.nin, 1)
        self.assertEqual(cf.nout, 1)
        self.assertEqual(cf.cname, 'f42_pypy_g_sqr')

    def test_genufunc(self):
        def foo(x):
            return x + 17
        uf = genufunc(foo, [
                (float, float),
                (int, int),
                ])
        self.assertEqual(uf(4), 21)
        x = array([1.1, 2.3])
        y = uf(x)
        self.assert_(allclose(y, [18.1, 19.3]))
        self.assert_(str(y.dtype).startswith('float'))
        
        x = array([1, 4])
        y = uf(x)
        self.assertEqual(list(y), [18, 21])
        self.assert_(str(y.dtype).startswith('int'))


class Arg_Tests(unittest.TestCase):
    
    def check_ufunc(self, f):
        for arg in (array([0.0, 1.0, 2.5]),
                    [0.0, 1.0, 2.5],
                    (0.0, 1.0, 2.5)):
            self.assert_(allclose(f(arg), [0.0, 1.0, 6.25]))
            
        self.assertEqual(f(3), 9)
        self.assert_(f(-2.5) - 6.25 < 1E-10)

    def test_direct(self):
        @mkufunc
        def f(x):
            return x * x
        self.check_ufunc(f)
        
    def test_noargs(self):
        @mkufunc()
        def f(x):
            return x * x
        self.check_ufunc(f)
        
    def test_varargs(self):
        for arg in (float,
                    [float],
                    [(float, float)]):
            @mkufunc(arg)
            def f(x):
                return x * x
            self.check_ufunc(f)


class Math_Tests(unittest.TestCase):
    
    def test_func1arg(self):
        for f in (math.exp, math.log, math.sqrt,
                  math.acos, math.asin, math.atan,
                  math.cos, math.sin, math.tan):
            @mkufunc
            def uf(x):
                return f(x)
            x = 0.4376
            a = uf(x)
            b = f(x)
            self.assert_(abs(a - b) < 1E-10, '%r %s != %s' % (f, a, b))
            xx = arange(0.1, 0.9, 0.01)
            a = uf(xx)
            b = [f(x) for x in xx]
            self.assert_(allclose(a, b))

    def test_arithmetic(self):
        def f(x):
            return (4 * x + 2) / (x * x - 7 * x + 1)
        uf = mkufunc(f)
        x = arange(0, 2, 0.1)
        self.assert_(allclose(uf(x), f(x)))
    

class Loop_Tests(unittest.TestCase):
    pass

class Switch_Tests(unittest.TestCase):
    pass

class FreeVariable_Tests(unittest.TestCase):
    pass

class Misc_Tests(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
