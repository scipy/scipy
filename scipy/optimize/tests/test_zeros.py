#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from math import sqrt, exp, sin, cos

from numpy.testing import TestCase, assert_almost_equal, assert_warns, \
                            assert_, run_module_suite, assert_allclose

from scipy.optimize import zeros as cc
from scipy.optimize import zeros

# Import testing parameters
from scipy.optimize._tstutils import functions, fstrings

class TestBasic(TestCase) :
    def run_check(self, method, name):
        a = .5
        b = sqrt(3)
        for function, fname in zip(functions, fstrings):
            zero, r = method(function, a, b, xtol=0.1e-12, full_output=True)
            assert_(r.converged)
            assert_almost_equal(zero, 1.0, decimal=12,
                err_msg='method %s, function %s' % (name, fname))

    def test_bisect(self):
        self.run_check(cc.bisect, 'bisect')
    def test_ridder(self):
        self.run_check(cc.ridder, 'ridder')
    def test_brentq(self):
        self.run_check(cc.brentq, 'brentq')
    def test_brenth(self):
        self.run_check(cc.brenth, 'brenth')

    def test_newton(self):
        f1 = lambda x: x**2 - 2*x - 1
        f1_1 = lambda x: 2*x - 2
        f1_2 = lambda x: 2.0 + 0*x

        f2 = lambda x: exp(x) - cos(x)
        f2_1 = lambda x: exp(x) + sin(x)
        f2_2 = lambda x: exp(x) + cos(x)

        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            x = zeros.newton(f, 3, tol=1e-6)
            assert_allclose(f(x), 0, atol=1e-6)
            x = zeros.newton(f, 3, fprime=f_1, tol=1e-6)
            assert_allclose(f(x), 0, atol=1e-6)
            x = zeros.newton(f, 3, fprime=f_1, fprime2=f_2, tol=1e-6)
            assert_allclose(f(x), 0, atol=1e-6)

    def test_deriv_zero_warning(self):
        func = lambda x: x**2
        dfunc = lambda x: 2*x
        assert_warns(RuntimeWarning, cc.newton, func, 0.0, dfunc)

if __name__ == '__main__' :
    run_module_suite()
