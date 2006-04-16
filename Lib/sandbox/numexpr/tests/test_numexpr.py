from numpy import *
from numpy.testing import *

set_package_path()
from numexpr import E, numexpr, evaluate
restore_path()

class test_numexpr(NumpyTestCase):
    def check_simple(self):
        ex = 2.0 * E.a + 3.0 * E.b * E.c
        func = numexpr(ex, signature=[('a', float), ('b', float), ('c', float)])
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

    def check_complex_expr(self):
        def complex(a, b, complex=__builtins__.complex):
            c = zeros(a.shape, dtype=complex)
            c.real = a
            c.imag = b
            return c
        a = arange(1e4)
        b = arange(1e4)**1e-5
        z = a + 1j*b
        x = z.imag
        x = sin(complex(a, b)).real + z.imag
        y = evaluate("sin(complex(a, b)).real + z.imag")
        assert_array_almost_equal(x, y)                                                                                                

tests = [
('MISC', ['b*c+d*e',
          '2*a+3*b',
          'sinh(a)',
          '2*a + (cos(3)+5)*sinh(cos(b))',
          '2*a + arctan2(a, b)',
          'where(a, 2, b)',
          'where((a-10).real, a, 2)',
          'cos(1+1)',
          '1+1',
          '1',
          'cos(a2)',
          '(a+1)**0'])]
optests = []
for op in list('+-*/%') + ['**']:
    optests.append("(a+1) %s (b+3)" % op)
    optests.append("3 %s (b+3)" % op)
    optests.append("(a+1) %s 4" % op)
    optests.append("2 %s (b+3)" % op)
    optests.append("(a+1) %s 2" % op)
    optests.append("(a+1) %s -1" % op)
    optests.append("(a+1) %s 0.5" % op)

tests.append(('OPERATIONS', optests))
cmptests = []
for op in ['<', '<=', '==', '>=', '>', '!=']:
    cmptests.append("a/2+5 %s b" % op)
    cmptests.append("a/2+5 %s 7" % op)
    cmptests.append("7 %s b" % op)
tests.append(('COMPARISONS', cmptests))
func1tests = []
for func in ['copy', 'ones_like', 'sin', 'cos', 'tan', 'sqrt', 'sinh', 'cosh', 'tanh']:
    func1tests.append("a + %s(b+c)" % func)
tests.append(('1-ARG FUNCS', func1tests))
func2tests = []
for func in ['arctan2', 'fmod']:
    func2tests.append("a + %s(b+c, d+1)" % func)
    func2tests.append("a + %s(b+c, 1)" % func)
    func2tests.append("a + %s(1, d+1)" % func)
tests.append(('2-ARG FUNCS', func2tests))
powtests = []
for n in (-2.5, -1.5, -1.3, -.5, 0, 0.5, 1, 0.5, 1, 2.3, 2.5):
    powtests.append("(a+1)**%s" % n)
tests.append(('POW TESTS', powtests))

def equal(a, b, exact):
    if exact:
        return (shape(a) == shape(b)) and alltrue(ravel(a) == ravel(b))
    else:
        return (shape(a) == shape(b)) and (allclose(ravel(a), ravel(b)) or alltrue(ravel(a) == ravel(b))) # XXX report a bug?

class Skip(Exception): pass

class test_expressions(NumpyTestCase):
    def check_expressions(self):
        for test_scalar in [0,1,2]:
            for dtype in [int, float, complex]:
                array_size = 100
                a = arange(array_size, dtype=dtype)
                a2 = zeros([array_size, array_size], dtype=dtype)
                b = arange(array_size, dtype=dtype) / array_size
                c = arange(array_size, dtype=dtype)
                d = arange(array_size, dtype=dtype)
                e = arange(array_size, dtype=dtype)
                if dtype == complex:
                    a = a.real
                    for x in [a2, b, c, d, e]:
                        x += 1j
                        x *= 1+1j
                if test_scalar == 1:
                    a = a[array_size/2] 
                if test_scalar == 2:
                    b = b[array_size/2]
                for optimization, exact in [('none', False), ('moderate', False), ('aggressive', False)]:
                    for section_name, section_tests in tests:
                        for expr in section_tests:
                            if dtype == complex and (
                                   '<' in expr or '>' in expr or '%' in expr
                                   or "arctan2" in expr or "fmod" in expr):
                                continue # skip complex comparisons
                            try:
                                try:
                                    npval = eval(expr)
                                except:
                                    raise Skip()
                                neval = evaluate(expr, optimization=optimization)
                                assert equal(npval, neval, exact), "%s (%s, %s, %s, %s)" % (expr, test_scalar, dtype.__name__, optimization, exact)
                            except Skip:
                                pass
                            except AssertionError:
                                raise
                            except NotImplementedError:
                                self.warn('%r not implemented for %s' % (expr,dtype.__name__))
                            except:
                                self.warn('numexpr error for expression %r' % (expr,))
                                raise

if __name__ == '__main__':
    NumpyTest().run()
