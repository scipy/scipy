#!/usr/bin/env python

from math import sqrt

from numpy.testing import TestCase, assert_almost_equal, assert_warns, \
                            assert_, run_module_suite

from scipy.optimize import zeros as cc

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

    def test_deriv_zero_warning(self):
        func = lambda x: x**2
        dfunc = lambda x: 2*x
        assert_warns(RuntimeWarning, cc.newton, func, 0.0, dfunc)

if __name__ == '__main__' :
    run_module_suite()
