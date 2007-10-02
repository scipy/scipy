import new
from numpy import *
from numpy.testing import *

set_package_path()
from numexpr import E, numexpr, evaluate, disassemble
restore_path()

class TestNumExpr(NumpyTestCase):
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

    def check_reductions(self):
        # Check that they compile OK.
        assert_equal(disassemble(numexpr("sum(x**2+2, axis=None)", [('x', float)])),
                    [('mul_fff', 't3', 'r1[x]', 'r1[x]'),
                     ('add_fff', 't3', 't3', 'c2[2.0]'),
                     ('sum_ffn', 'r0', 't3', None)])
        assert_equal(disassemble(numexpr("sum(x**2+2, axis=1)", [('x', float)])),
                    [('mul_fff', 't3', 'r1[x]', 'r1[x]'),
                     ('add_fff', 't3', 't3', 'c2[2.0]'),
                     ('sum_ffn', 'r0', 't3', 1)])
        assert_equal(disassemble(numexpr("prod(x**2+2, axis=2)", [('x', float)])),
                    [('mul_fff', 't3', 'r1[x]', 'r1[x]'),
                     ('add_fff', 't3', 't3', 'c2[2.0]'),
                     ('prod_ffn', 'r0', 't3', 2)])
        # Check that full reductions work.
        x = arange(10.0)
        assert_equal(evaluate("sum(x**2+2,axis=0)"), sum(x**2+2,axis=0))
        assert_equal(evaluate("prod(x**2+2,axis=0)"), prod(x**2+2,axis=0))
        # Check that reductions along an axis work
        y = arange(9.0).reshape(3,3)
        assert_equal(evaluate("sum(y**2, axis=1)"), sum(y**2, axis=1))
        assert_equal(evaluate("sum(y**2, axis=0)"), sum(y**2, axis=0))
        assert_equal(evaluate("sum(y**2, axis=None)"), sum(y**2, axis=None))
        assert_equal(evaluate("prod(y**2, axis=1)"), prod(y**2, axis=1))
        assert_equal(evaluate("prod(y**2, axis=0)"), prod(y**2, axis=0))
        assert_equal(evaluate("prod(y**2, axis=None)"), prod(y**2, axis=None))
        # Check integers
        x = x.astype(int)
        assert_equal(evaluate("sum(x**2+2,axis=0)"), sum(x**2+2,axis=0))
        assert_equal(evaluate("prod(x**2+2,axis=0)"), prod(x**2+2,axis=0))
        # Check complex
        x = x + 5j
        assert_equal(evaluate("sum(x**2+2,axis=0)"), sum(x**2+2,axis=0))
        assert_equal(evaluate("prod(x**2+2,axis=0)"), prod(x**2+2,axis=0))
        # Check boolean (should cast to integer)
        x = (arange(10) % 2).astype(bool)
        assert_equal(evaluate("prod(x,axis=0)"), prod(x,axis=0))
        assert_equal(evaluate("sum(x,axis=0)"), sum(x,axis=0))

    def check_axis(self):
        y = arange(9.0).reshape(3,3)
        try:
            evaluate("sum(y, axis=2)")
        except ValueError:
            pass
        else:
            raise ValueError("should raise exception!")
        try:
            evaluate("sum(y, axis=-3)")
        except ValueError:
            pass
        else:
            raise ValueError("should raise exception!")




    def check_r0_reuse(self):
        assert_equal(disassemble(numexpr("x**2+2", [('x', float)])),
                    [('mul_fff', 'r0', 'r1[x]', 'r1[x]'),
                     ('add_fff', 'r0', 'r0', 'c2[2.0]')])

class TestEvaluate(NumpyTestCase):
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
        def complex(a, b):
            c = zeros(a.shape, dtype=complex_)
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


    def check_complex_strides(self):
        a = arange(100).reshape(10,10)[::2]
        b = arange(50).reshape(5,10)
        assert_array_equal(evaluate("a+b"), a+b)
        c = empty([10], dtype=[('c1', int32), ('c2', uint16)])
        c['c1'] = arange(10)
        c['c2'].fill(0xaaaa)
        c1 = c['c1']
        a0 = a[0]
        assert_array_equal(evaluate("c1"), c1)
        assert_array_equal(evaluate("a0+c1"), a0+c1)


    def check_broadcasting(self):
        a = arange(100).reshape(10,10)[::2]
        c = arange(10)
        d = arange(5).reshape(5,1)
        assert_array_equal(evaluate("a+c"), a+c)
        assert_array_equal(evaluate("a+d"), a+d)
        expr = numexpr("2.0*a+3.0*c",[('a',float),('c', float)])
        assert_array_equal(expr(a,c), 2.0*a+3.0*c)

    def check_all_scalar(self):
        a = 3.
        b = 4.
        assert_equal(evaluate("a+b"), a+b)
        expr = numexpr("2*a+3*b",[('a',float),('b', float)])
        assert_equal(expr(a,b), 2*a+3*b)

    def check_run(self):
        a = arange(100).reshape(10,10)[::2]
        b = arange(10)
        expr = numexpr("2*a+3*b",[('a',float),('b', float)])
        assert_array_equal(expr(a,b), expr.run(a,b))

    def check_illegal_value(self):
        a = arange(3)
        try:
            evaluate("a < [0, 0, 0]")
        except TypeError:
            pass
        else:
            self.fail()


tests = [
('MISC', ['b*c+d*e',
          '2*a+3*b',
          'sinh(a)',
          '2*a + (cos(3)+5)*sinh(cos(b))',
          '2*a + arctan2(a, b)',
          'arcsin(0.5)',
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
    cmptests.append("7.0 %s 5" % op)
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
        return (shape(a) == shape(b)) and alltrue(ravel(a) == ravel(b),axis=0)
    else:
        return (shape(a) == shape(b)) and (allclose(ravel(a), ravel(b)) or alltrue(ravel(a) == ravel(b),axis=0)) # XXX report a bug?

class Skip(Exception): pass

class TestExpressions(NumpyTestCase):
    pass

def generate_check_expressions():
    test_no = [0]
    def make_check_method(a, a2, b, c, d, e, x, expr,
                          test_scalar, dtype, optimization, exact):
        this_locals = locals()
        def method(self):
            try:
                npval = eval(expr, globals(), this_locals)
            except:
                return
            try:
                neval = evaluate(expr, local_dict=this_locals,
                                 optimization=optimization)
                assert equal(npval, neval, exact), \
                    """%r
(test_scalar=%r, dtype=%r, optimization=%r, exact=%r,
 npval=%r (%r), neval=%r (%r))""" % (expr, test_scalar, dtype.__name__,
                                     optimization, exact,
                                     npval, type(npval), neval, type(neval))
            except AssertionError:
                raise
            except NotImplementedError:
                self.warn('%r not implemented for %s' % (expr,dtype.__name__))
            except:
                self.warn('numexpr error for expression %r' % (expr,))
                raise
        test_no[0] += 1
        name = 'check_%04d' % (test_no[0],)
        setattr(test_expressions, name,
                new.instancemethod(method, None, test_expressions))
    x = None
    for test_scalar in [0,1,2]:
        for dtype in [int, float, complex]:
            array_size = 100
            a = arange(2*array_size, dtype=dtype)[::2]
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
                        if dtype == int and test_scalar and expr == '(a+1) ** -1':
                            continue
                        make_check_method(a, a2, b, c, d, e, x,
                                          expr, test_scalar, dtype,
                                          optimization, exact)

generate_check_expressions()

if __name__ == '__main__':
    NumpyTest().run()
