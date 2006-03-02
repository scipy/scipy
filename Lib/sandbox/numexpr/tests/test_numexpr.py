from numpy import array, arange
from numpy.testing import *

set_package_path()
from numexpr import E, numexpr, evaluate
restore_path()

class test_numexpr(NumpyTestCase):
    def check_simple(self):
        ex = 2.0 * E.a + 3.0 * E.b * E.c
        func = numexpr(ex, input_order=('a', 'b', 'c'))
        x = func(array([1., 2, 3]), array([4., 5, 6]), array([7., 8, 9]))
        assert_array_equal(x, array([  86.,  124.,  168.]))

    def check_simple_expr_small_array(self):
        func = numexpr(E.a)
        x = arange(100.0)
        y = func(x)
        assert_array_equal(x, y)

    def check_simple_expr(self):
        func = numexpr(E.a)
        x = arange(1e5)
        y = func(x)
        assert_array_equal(x, y)

    def check_rational_expr(self):
        func = numexpr((E.a + 2.0*E.b) / (1 + E.a + 4*E.b*E.b))
        a = arange(1e5)
        b = arange(1e5) * 0.1
        x = (a + 2*b) / (1 + a + 4*b*b)
        y = func(a, b)
        assert_array_equal(x, y)

class test_evaluate(NumpyTestCase):
    def check_simple(self):
        a = array([1., 2., 3.])
        b = array([4., 5., 6.])
        c = array([7., 8., 9.])
        x = evaluate("2*a + 3*b*c")
        assert_array_equal(x, array([  86.,  124.,  168.]))

    def check_simple_expr_small_array(self):
        x = arange(100.0)
        y = evaluate("x")
        assert_array_equal(x, y)

    def check_simple_expr(self):
        x = arange(1e5)
        y = evaluate("x")
        assert_array_equal(x, y)

    def check_rational_expr(self):
        a = arange(1e5)
        b = arange(1e5) * 0.1
        x = (a + 2*b) / (1 + a + 4*b*b)
        y = evaluate("(a + 2*b) / (1 + a + 4*b*b)")
        assert_array_equal(x, y)

if __name__ == '__main__':
    NumpyTest().run()
