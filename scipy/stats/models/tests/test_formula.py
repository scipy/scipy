"""
Test functions for models.formula
"""

import string

import numpy as N
import numpy.random as R
import numpy.linalg as L
from scipy.testing import *

from scipy.stats.models import utils, formula, contrast

class TestTerm(TestCase):

    def test_init(self):
        t1 = formula.Term("trivial")
        sqr = lambda x: x*x

        t2 = formula.Term("not_so_trivial", sqr, "sqr")

        self.assertRaises(ValueError, formula.Term, "name", termname=0)


    def test_str(self):
        t = formula.Term("name")
        s = str(t)

    def test_add(self):
        t1 = formula.Term("t1")
        t2 = formula.Term("t2")
        f = t1 + t2
        self.assert_(isinstance(f, formula.Formula))
        self.assert_(f.hasterm(t1))
        self.assert_(f.hasterm(t2))

    def test_mul(self):
        t1 = formula.Term("t1")
        t2 = formula.Term("t2")
        f = t1 * t2
        self.assert_(isinstance(f, formula.Formula))

        intercept = formula.Term("intercept")
        f = t1 * intercept
        self.assertEqual(str(f), str(formula.Formula(t1)))

        f = intercept * t1
        self.assertEqual(str(f), str(formula.Formula(t1)))

class TestFormula(TestCase):

    def setUp(self):
        self.X = R.standard_normal((40,10))
        self.namespace = {}
        self.terms = []
        for i in range(10):
            name = '%s' % string.uppercase[i]
            self.namespace[name] = self.X[:,i]
            self.terms.append(formula.Term(name))

        self.formula = self.terms[0]
        for i in range(1, 10):
            self.formula += self.terms[i]
        self.formula.namespace = self.namespace

    def test_namespace(self):
        space1 = {'X':N.arange(50), 'Y':N.arange(50)*2}
        space2 = {'X':N.arange(20), 'Y':N.arange(20)*2}
        X = formula.Term('X')
        Y = formula.Term('Y')

        X.namespace = space1
        assert_almost_equal(X(), N.arange(50))

        Y.namespace = space2
        assert_almost_equal(Y(), N.arange(20)*2)

        f = X + Y

        f.namespace = space1
        self.assertEqual(f().shape, (2,50))
        assert_almost_equal(Y(), N.arange(50)*2)
        assert_almost_equal(X(), N.arange(50))

        f.namespace = space2
        self.assertEqual(f().shape, (2,20))
        assert_almost_equal(Y(), N.arange(20)*2)
        assert_almost_equal(X(), N.arange(20))


    def test_termcolumns(self):
        t1 = formula.Term("A")
        t2 = formula.Term("B")
        f = t1 + t2 + t1 * t2
        def other(val):
            return N.array([3.2*val,4.342*val**2, 5.234*val**3])
        q = formula.Quantitative(['other%d' % i for i in range(1,4)], termname='other', func=t1, transform=other)
        f += q
        q.namespace = f.namespace = self.formula.namespace
        assert_almost_equal(q(), f()[f.termcolumns(q)])


    def test_str(self):
        s = str(self.formula)

    def test_call(self):
        x = self.formula()
        self.assertEquals(N.array(x).shape, (10, 40))

    def test_design(self):
        x = self.formula.design()
        self.assertEquals(x.shape, (40, 10))

    def test_product(self):
        prod = self.terms[0] * self.terms[2]
        self.formula += prod
        x = self.formula.design()
        p = self.formula['A*C']
        col = self.formula.termcolumns(prod, dict=False)
        assert_almost_equal(N.squeeze(x[:,col]), self.X[:,0] * self.X[:,2])
        assert_almost_equal(N.squeeze(p()), self.X[:,0] * self.X[:,2])

    def test_intercept1(self):
        prod = self.terms[0] * self.terms[2]
        self.formula += formula.I
        icol = self.formula.names().index('intercept')
        assert_almost_equal(self.formula()[icol], N.ones((40,)))

    def test_intercept2(self):
        prod = self.terms[0] * self.terms[2]
        self.formula += formula.I
        icol = self.formula.names().index('intercept')
        assert_almost_equal(self.formula()[icol], N.ones((40,)))

    def test_intercept3(self):
        prod = self.terms[0] * formula.I
        prod.namespace = self.formula.namespace
        assert_almost_equal(N.squeeze(prod()), self.terms[0]())



    def test_contrast1(self):
        term = self.terms[0] + self.terms[2]
        c = contrast.Contrast(term, self.formula)
        c.getmatrix()
        col1 = self.formula.termcolumns(self.terms[0], dict=False)
        col2 = self.formula.termcolumns(self.terms[1], dict=False)
        test = [[1] + [0]*9, [0]*2 + [1] + [0]*7]
        assert_almost_equal(c.matrix, test)

    def test_contrast2(self):

        dummy = formula.Term('zero')
        self.namespace['zero'] = N.zeros((40,), N.float64)
        term = dummy + self.terms[2]
        c = contrast.Contrast(term, self.formula)
        c.getmatrix()
        test = [0]*2 + [1] + [0]*7
        assert_almost_equal(c.matrix, test)

    def test_contrast3(self):

        X = self.formula.design()
        P = N.dot(X, L.pinv(X))

        dummy = formula.Term('noise')
        resid = N.identity(40) - P
        self.namespace['noise'] = N.transpose(N.dot(resid, R.standard_normal((40,5))))
        terms = dummy + self.terms[2]
        terms.namespace = self.formula.namespace
        c = contrast.Contrast(terms, self.formula)
        c.getmatrix()
        self.assertEquals(c.matrix.shape, (10,))

    def test_power(self):

        t = self.terms[2]
        t2 = t**2
        t.namespace = t2.namespace = self.formula.namespace
        assert_almost_equal(t()**2, t2())

    def test_quantitative(self):
        t = self.terms[2]
        sint = formula.Quantitative('t', func=t, transform=N.sin)
        t.namespace = sint.namespace = self.formula.namespace
        assert_almost_equal(N.sin(t()), sint())

    def test_factor1(self):
        f = ['a','b','c']*10
        fac = formula.Factor('ff', f)
        fac.namespace = {'ff':f}
        self.assertEquals(list(fac.values()), f)

    def test_factor2(self):
        f = ['a','b','c']*10
        fac = formula.Factor('ff', f)
        fac.namespace = {'ff':f}
        self.assertEquals(fac().shape, (3,30))

    def test_factor3(self):
        f = ['a','b','c']*10
        fac = formula.Factor('ff', f)
        fac.namespace = {'ff':f}
        m = fac.main_effect(reference=1)
        self.assertEquals(m().shape, (2,30))

    def test_factor4(self):
        f = ['a','b','c']*10
        fac = formula.Factor('ff', f)
        fac.namespace = {'ff':f}
        m = fac.main_effect(reference=2)
        r = N.array([N.identity(3)]*10)
        r.shape = (30,3)
        r = r.T
        _m = N.array([r[0]-r[2],r[1]-r[2]])
        assert_almost_equal(_m, m())

    def test_contrast4(self):

        f = self.formula + self.terms[5] + self.terms[5]

        estimable = False

        c = contrast.Contrast(self.terms[5], f)
        c.getmatrix()

        self.assertEquals(estimable, False)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
