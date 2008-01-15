import new
from numpy import array, arange, sin, zeros, sum, int32, empty, \
     prod, uint16, complex_, float64, rec
from scipy.testing import *


from scipy.sandbox.numexpr import E, numexpr, evaluate, disassemble


class test_numexpr(TestCase):
    def test_simple(self):
        ex = 2.0 * E.a + 3.0 * E.b * E.c
        func = numexpr(ex, signature=[('a', float), ('b', float), ('c', float)])
        x = func(array([1., 2, 3]), array([4., 5, 6]), array([7., 8, 9]))
        assert_array_equal(x, array([  86.,  124.,  168.]))

    def test_simple_expr_small_array(self):
        func = numexpr(E.a)
        x = arange(100.0)
        y = func(x)
        assert_array_equal(x, y)

    def test_simple_expr(self):
        func = numexpr(E.a)
        x = arange(1e5)
        y = func(x)
        assert_array_equal(x, y)

    def test_rational_expr(self):
        func = numexpr((E.a + 2.0*E.b) / (1 + E.a + 4*E.b*E.b))
        a = arange(1e5)
        b = arange(1e5) * 0.1
        x = (a + 2*b) / (1 + a + 4*b*b)
        y = func(a, b)
        assert_array_equal(x, y)

    def test_reductions(self):
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
        # Check longs
        x = x.astype(long)
        assert_equal(evaluate("sum(x**2+2,axis=0)"), sum(x**2+2,axis=0))
        assert_equal(evaluate("prod(x**2+2,axis=0)"), prod(x**2+2,axis=0))
        # Check complex
        x = x + 5j
        assert_equal(evaluate("sum(x**2+2,axis=0)"), sum(x**2+2,axis=0))
        assert_equal(evaluate("prod(x**2+2,axis=0)"), prod(x**2+2,axis=0))

    def test_axis(self):
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




    def test_r0_reuse(self):
        assert_equal(disassemble(numexpr("x**2+2", [('x', float)])),
                    [('mul_fff', 'r0', 'r1[x]', 'r1[x]'),
                     ('add_fff', 'r0', 'r0', 'c2[2.0]')])

class test_evaluate(TestCase):
    def test_simple(self):
        a = array([1., 2., 3.])
        b = array([4., 5., 6.])
        c = array([7., 8., 9.])
        x = evaluate("2*a + 3*b*c")
        assert_array_equal(x, array([  86.,  124.,  168.]))

    def test_simple_expr_small_array(self):
        x = arange(100.0)
        y = evaluate("x")
        assert_array_equal(x, y)

    def test_simple_expr(self):
        x = arange(1e5)
        y = evaluate("x")
        assert_array_equal(x, y)

    def test_rational_expr(self):
        a = arange(1e5)
        b = arange(1e5) * 0.1
        x = (a + 2*b) / (1 + a + 4*b*b)
        y = evaluate("(a + 2*b) / (1 + a + 4*b*b)")
        assert_array_equal(x, y)

    def test_complex_expr(self):
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


    def test_complex_strides(self):
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


    def test_broadcasting(self):
        a = arange(100).reshape(10,10)[::2]
        c = arange(10)
        d = arange(5).reshape(5,1)
        assert_array_equal(evaluate("a+c"), a+c)
        assert_array_equal(evaluate("a+d"), a+d)
        expr = numexpr("2.0*a+3.0*c",[('a',float),('c', float)])
        assert_array_equal(expr(a,c), 2.0*a+3.0*c)

    def test_all_scalar(self):
        a = 3.
        b = 4.
        assert_equal(evaluate("a+b"), a+b)
        expr = numexpr("2*a+3*b",[('a',float),('b', float)])
        assert_equal(expr(a,b), 2*a+3*b)

    def test_run(self):
        a = arange(100).reshape(10,10)[::2]
        b = arange(10)
        expr = numexpr("2*a+3*b",[('a',float),('b', float)])
        assert_array_equal(expr(a,b), expr.run(a,b))

    def test_illegal_value(self):
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
          'where(a != 0.0, 2, b)',
          'where((a-10).real != 0.0, a, 2)',
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

class test_expressions(TestCase):
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
        for dtype in [int, long, float, complex]:
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
                        if dtype in (int, long) and test_scalar and expr == '(a+1) ** -1':
                            continue
                        make_check_method(a, a2, b, c, d, e, x,
                                          expr, test_scalar, dtype,
                                          optimization, exact)

generate_check_expressions()

class test_int32_int64(TestCase):
    def test_small_long(self):
        # Small longs should not be downgraded to ints.
        res = evaluate('42L')
        assert_array_equal(res, 42)
        self.assertEqual(res.dtype.name, 'int64')

    def test_big_int(self):
        # Big ints should be promoted to longs.
        # This test may only fail under 64-bit platforms.
        res = evaluate('2**40')
        assert_array_equal(res, 2**40)
        self.assertEqual(res.dtype.name, 'int64')

    def test_long_constant_promotion(self):
        int32array = arange(100, dtype='int32')
        res = int32array * 2
        res32 = evaluate('int32array * 2')
        res64 = evaluate('int32array * 2L')
        assert_array_equal(res, res32)
        assert_array_equal(res, res64)
        self.assertEqual(res32.dtype.name, 'int32')
        self.assertEqual(res64.dtype.name, 'int64')

    def test_int64_array_promotion(self):
        int32array = arange(100, dtype='int32')
        int64array = arange(100, dtype='int64')
        respy = int32array * int64array
        resnx = evaluate('int32array * int64array')
        assert_array_equal(respy, resnx)
        self.assertEqual(resnx.dtype.name, 'int64')

class test_strings(TestCase):
    BLOCK_SIZE1 = 128
    BLOCK_SIZE2 = 8
    str_list1 = ['foo', 'bar', '', '  ']
    str_list2 = ['foo', '', 'x', ' ']
    str_nloops = len(str_list1) * (BLOCK_SIZE1 + BLOCK_SIZE2 + 1)
    str_array1 = array(str_list1 * str_nloops)
    str_array2 = array(str_list2 * str_nloops)
    str_constant = 'doodoo'

    def test_null_chars(self):
        str_list = [
            '\0\0\0', '\0\0foo\0', '\0\0foo\0b', '\0\0foo\0b\0',
            'foo\0', 'foo\0b', 'foo\0b\0', 'foo\0bar\0baz\0\0' ]
        for s in str_list:
            r = evaluate('s')
            self.assertEqual(s, r.tostring())  # check *all* stored data

    def test_compare_copy(self):
        sarr = self.str_array1
        expr = 'sarr'
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_compare_array(self):
        sarr1 = self.str_array1
        sarr2 = self.str_array2
        expr = 'sarr1 >= sarr2'
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_compare_variable(self):
        sarr = self.str_array1
        svar = self.str_constant
        expr = 'sarr >= svar'
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_compare_constant(self):
        sarr = self.str_array1
        expr = 'sarr >= %r' % self.str_constant
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_add_string_array(self):
        sarr1 = self.str_array1
        sarr2 = self.str_array2
        expr = 'sarr1 + sarr2'
        self.assert_missing_op('add_sss', expr, locals())

    def test_add_numeric_array(self):
        sarr = self.str_array1
        narr = arange(len(sarr), dtype='int32')
        expr = 'sarr >= narr'
        self.assert_missing_op('ge_bsi', expr, locals())

    def assert_missing_op(self, op, expr, local_dict):
        msg = "expected NotImplementedError regarding '%s'" % op
        try:
            evaluate(expr, local_dict)
        except NotImplementedError, nie:
            if "'%s'" % op not in nie.args[0]:
                self.fail(msg)
        else:
            self.fail(msg)

    def test_compare_prefix(self):
        # Check comparing two strings where one is a prefix of the
        # other.
        for s1, s2 in [ ('foo', 'foobar'), ('foo', 'foo\0bar'),
                        ('foo\0a', 'foo\0bar') ]:
            self.assert_(evaluate('s1 < s2'))
            self.assert_(evaluate('s1 <= s2'))
            self.assert_(evaluate('~(s1 == s2)'))
            self.assert_(evaluate('~(s1 >= s2)'))
            self.assert_(evaluate('~(s1 > s2)'))

        # Check for NumPy array-style semantics in string equality.
        s1, s2 = 'foo', 'foo\0\0'
        self.assert_(evaluate('s1 == s2'))

# Case for testing selections in fields which are aligned but whose
# data length is not an exact multiple of the length of the record.
# The following test exposes the problem only in 32-bit machines,
# because in 64-bit machines 'c2' is unaligned.  However, this should
# check most platforms where, while not unaligned, 'len(datatype) >
# boundary_alignment' is fullfilled.
class test_irregular_stride(TestCase):
    def test_select(self):
        f0 = arange(10, dtype=int32)
        f1 = arange(10, dtype=float64)

        irregular = rec.fromarrays([f0, f1])

        f0 = irregular['f0']
        f1 = irregular['f1']

        i0 = evaluate('f0 < 5')
        i1 = evaluate('f1 < 5')

        assert_array_equal(f0[i0], arange(5, dtype=int32))
        assert_array_equal(f1[i1], arange(5, dtype=float64))


if __name__ == '__main__':
    nose.run(argv=['', __file__])
