from numpy import *
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


tests = [
('MISC', ['b*c+d*e',
          '2*a+3*b',
          'sinh(a)',
          '2*a + (cos(3)+5)*sinh(cos(b))',
          '2*a + arctan2(a, b)',
          'where(0.1*a > arctan2(a, b), 2*a, arctan2(a,b))',
          'where(a, 2, b)',
          'where(a-10, a, 2)',
          'cos(1+1)',
          '1+1',
          '1',
          'cos(a2)'])]
optests = []
for op in list('+-*/%') + ['**']:
    optests.append("(a+1) %s (b+3)" % op)
    optests.append("3 %s (b+3)" % op)
    optests.append("(a+1) %s 4" % op)
tests.append(('OPERATIONS', optests))
cmptests = []
for op in ['<', '<=', '==', '>=', '>', '!=']:
    cmptests.append("a/2+5 %s b" % op)
tests.append(('COMPARISONS', cmptests))
func1tests = []
for func in ['sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh']:
    func1tests.append("a + %s(b+c)" % func)
tests.append(('1-ARG FUNCS', func1tests))
func2tests = []
for func in ['arctan2', 'fmod']:
    func2tests.append("a + %s(b+c, d+1)" % func)
    func2tests.append("a + %s(b+c, 1)" % func)
    func2tests.append("a + %s(1, d+1)" % func)
tests.append(('2-ARG FUNCS', func2tests))

class test_expressions(NumpyTestCase):
    def check_expressions(self):
        array_size = 1e2
        a = arange(array_size)
        a2 = zeros([array_size, array_size])
        b = arange(array_size)
        c = arange(array_size)
        d = arange(array_size)
        e = arange(array_size)

        try:
            for section_name, section_tests in tests:
                for expr in section_tests:
                    npval = eval(expr)
                    neval = evaluate(expr)
                    assert shape(npval) == shape(neval), expr
                    assert alltrue(ravel(npval) == ravel(neval)), expr
        except AssertionError:
            raise
        except:
            self.warn('numexpr error for expression %r' % (expr,))
            raise

if __name__ == '__main__':
    NumpyTest().run()
