import math
import unittest

from numpy import array, arange, allclose

from api import Cfunc, genufunc, mkufunc


class Util:

    def assertClose(self, x, y):
        self.assert_(allclose(x, y), '%s != %s' % (x, y))
            

class Internal_Tests(unittest.TestCase, Util):
    
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
        self.assertClose(y, [18.1, 19.3])
        self.assert_(str(y.dtype).startswith('float'))
        
        x = array([1, 4])
        y = uf(x)
        self.assertEqual(list(y), [18, 21])
        self.assert_(str(y.dtype).startswith('int'))


class Arg_Tests(unittest.TestCase, Util):
    
    def check_ufunc(self, f):
        for arg in (array([0.0, 1.0, 2.5]),
                    [0.0, 1.0, 2.5],
                    (0.0, 1.0, 2.5)):
            self.assertClose(f(arg), [0.0, 1.0, 6.25])
            
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

    def test_int(self):
        @mkufunc(int)
        def f(x):
            return x * x
        self.assertEqual(f(3), 9)
        self.assert_(isinstance(f(42), int))
        
    def test_mixed(self):
        @mkufunc([(int, float, int), float])
        def f(n, x):
            return n + x * x
        
        y = f(2, 3.9)            # Note that int(2 + 3.9 * 3.9) = 17
        self.assertEqual(y, 17)
        self.assert_(isinstance(y, int))
        
        y = f(2.0, 3.9)
        self.assertClose(y, 17.21)
        self.assert_(isinstance(y, float))
        
    def test_exceptions(self):
        def f(x):
            return x

        self.assertRaises(TypeError, mkufunc, {})
        self.assertRaises(TypeError, mkufunc([(float,)]), f)
        self.assertRaises(TypeError, mkufunc([3*(float,)]), f)
        self.assertRaises(TypeError, mkufunc([{}]), f)
        self.assertRaises(TypeError, mkufunc([(int, {})]), f)
        self.assertRaises(ValueError, mkufunc([]), f)
        

class Math_Tests(unittest.TestCase, Util):
    
    def assertFuncsEqual(self, uf, f):
        x = 0.4376
        a = uf(x)
        b = f(x)
        self.assertClose(a, b)
        xx = arange(0.1, 0.9, 0.01)
        a = uf(xx)
        b = [f(x) for x in xx]
        self.assertClose(a, b)
        
    def test_pi(self):
        @mkufunc
        def f(x): return math.pi * x
        self.assertFuncsEqual(f, lambda x: math.pi * x)

    def test_e(self):
        @mkufunc#(show=1)
        def f(x): return math.e * x
        self.assertFuncsEqual(f, lambda x: math.e * x)
   
    def test_exp(self):
        @mkufunc
        def f(x): return math.exp(x)
        self.assertFuncsEqual(f, math.exp)

    def test_log(self):
        @mkufunc
        def f(x): return math.log(x)
        self.assertFuncsEqual(f, math.log)

    def test_log10(self):
        @mkufunc
        def f(x): return math.log10(x)
        self.assertFuncsEqual(f, math.log10)

    def test_sqrt(self):
        @mkufunc
        def f(x): return math.sqrt(x)
        self.assertFuncsEqual(f, math.sqrt)
        
    def test_cos(self):
        @mkufunc
        def f(x): return math.cos(x)
        self.assertFuncsEqual(f, math.cos)
        
    def test_sin(self):
        @mkufunc
        def f(x): return math.sin(x)
        self.assertFuncsEqual(f, math.sin)
        
    def test_tan(self):
        @mkufunc
        def f(x): return math.tan(x)
        self.assertFuncsEqual(f, math.tan)
        
    def test_cosh(self):
        @mkufunc
        def f(x): return math.cosh(x)
        self.assertFuncsEqual(f, math.cosh)

    def test_sinh(self):
        @mkufunc
        def f(x): return math.sinh(x)
        self.assertFuncsEqual(f, math.sinh)
        
    def test_tanh(self):
        @mkufunc
        def f(x): return math.tanh(x)
        self.assertFuncsEqual(f, math.tanh)
        
    def test_acos(self):
        @mkufunc
        def f(x): return math.acos(x)
        self.assertFuncsEqual(f, math.acos)

    def test_asin(self):
        @mkufunc
        def f(x): return math.asin(x)
        self.assertFuncsEqual(f, math.asin)

    def test_atan(self):
        @mkufunc
        def f(x): return math.atan(x)
        self.assertFuncsEqual(f, math.atan)

    def test_atan2(self):
        @mkufunc
        def f(x, y):
            return math.atan2(x, y)

        self.assertClose(f(4, 5), math.atan2(4, 5))
        
        xx = array([1.0, 3.0, -2.4,  3.1, -2.3, -1.0])
        yy = array([1.0, 2.0,  7.5, -8.7,  0.0, -3.2])
        a = f(xx, yy)
        b = [math.atan2(x, y) for x, y in zip(xx, yy)]
        self.assertClose(a, b)

    def test_pow(self):
        @mkufunc
        def f(x, y):
            return math.pow(x, y)
        
        self.assertClose(f(2, 3), 8)
        
        xx = array([1.0, 3.0, 2.4, 0.0, 0.0, 0.0, 2.3,  2.0])
        yy = array([1.0, 2.0, 7.5, 0.0, 0.5, 1.0, 0.0, -1.0])
        a = f(xx, yy)
        b = [math.pow(x, y) for x, y in zip(xx, yy)]
        self.assertClose(a, b)
       
    def test_arithmetic(self):
        def f(x):
            return (4 * x + 2) / (x * x - 7 * x + 1)
        uf = mkufunc(f)
        x = arange(0, 2, 0.1)
        self.assertClose(uf(x), f(x))


class Control_Flow_Tests(unittest.TestCase):

    def test_if(self):
        @mkufunc(int)
        def f(n):
            if n < 4:
                return n
            else:
                return n * n

        self.assertEqual(f(3), 3)
        self.assertEqual(f(4), 16)

    def test_switch(self):
        @mkufunc(int)
        def f(n):
            if n < 4:
                return n
            elif n == 4:
                return 42
            elif n == 5:
                return 73
            else:
                return n * n

        self.assertEqual(f(3), 3)
        self.assertEqual(f(4), 42)
        self.assertEqual(f(5), 73)
        self.assertEqual(f(6), 36)

    def test_loop(self):
        @mkufunc(int)
        def f(n):
            res = 0
            for i in xrange(n):
                res += i*i
            return res

        self.assertEqual(f(3), 5)
        self.assertEqual(f(95), 281295)


class FreeVariable_Tests(unittest.TestCase, Util):

    def test_const(self):
        a = 13.6
        @mkufunc
        def f(x):
            return a * x
        
        x = arange(0, 1, 0.1)
        self.assertClose(f(x), a * x)


class Misc_Tests(unittest.TestCase, Util):

    def test_lambda(self):
        self.assertRaises(AssertionError, mkufunc, lambda x: x*x + 2)
    
    def test_angle(self):
        from math import sin, pi, sqrt
        @mkufunc
        def sin_deg(angle):
            return sin(angle / 180.0 * pi)
        
        self.assertClose(sin_deg([0, 30, 45, 60, 90, 180, 270, 360]),
                         [0, 0.5, 1/sqrt(2), sqrt(3)/2, 1, 0, -1, 0])
  



if __name__ == '__main__':
    unittest.main()
