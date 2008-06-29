import math
import unittest

from numpy import array, arange, allclose

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
        self.assert_(abs(y - 17.21) < 1E-10)
        self.assert_(isinstance(y, float))
        
    def test_exceptions(self):
        def f(x):
            return x

        self.assertRaises(TypeError, mkufunc, {})
        self.assertRaises(TypeError, mkufunc([(float,)]), f)
        self.assertRaises(TypeError, mkufunc([3*(float,)]), f)
        self.assertRaises(TypeError, mkufunc([{}]), f)
        self.assertRaises(TypeError, mkufunc([(int, {})]), f)
        

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
            
    def test_func2arg(self):
        @mkufunc
        def f(x, y):
            return math.atan2(x, y)
        
        xx = array([1.0, 3.0, -2.4,  3.1, -2.3])
        yy = array([1.0, 2.0,  7.5, -8.7,  0.0])
        a = f(xx, yy)
        b = [math.atan2(x, y) for x, y in zip(xx, yy)]
        self.assert_(allclose(a, b))
        
    def test_arithmetic(self):
        def f(x):
            return (4 * x + 2) / (x * x - 7 * x + 1)
        uf = mkufunc(f)
        x = arange(0, 2, 0.1)
        self.assert_(allclose(uf(x), f(x)))
    
        def f(x, y, z):
            return x * y * z
        uf = mkufunc(f)
        x = arange(0, 1, 0.1)
        y = 2 * x
        z = 3 * x
        self.assert_(allclose(uf(x, y, z), f(x, y, z)))


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


class FreeVariable_Tests(unittest.TestCase):

    def test_const(self):
        a = 13.6
        @mkufunc
        def f(x):
            return a * x
        
        x = arange(0, 1, 0.1)
        self.assert_(allclose(f(x), a * x))

    def test_const2(self):
        from math import sin, pi, sqrt
        @mkufunc
        def sin_deg(angle):
            return sin(angle / 180.0 * pi)
        
        self.assert_(allclose(sin_deg([0, 30, 45, 60, 90, 180, 270, 360]),
                              [0, 0.5, 1/sqrt(2), sqrt(3)/2, 1, 0, -1, 0]))
        

class Misc_Tests(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
