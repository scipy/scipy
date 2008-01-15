#!/usr/bin/env python


from scipy.testing import *

from scipy.optimize import zeros as cc

from math import sin,sqrt,log
from random import random

def f1(x) :
    return x*(x-1.)

def f2(x) :
    return x**2 - 1

def f3(x) :
    return x*(x-1.)*(x-2.)*(x-3.)

def f4(x) :
    if x > 1 : return  1.0 + .1*x
    if x < 1 : return -1.0 + .1*x
    return 0

def f5(x) :
    if x != 1 : return 1.0/(1. - x)
    return 0

def f6(x) :
    if   x > 1 : return random()
    elif x < 1 : return -random()
    else : return 0

description = """
f2 is a symmetric parabola, x**2 - 1
f3 is a quartic polynomial with large hump in interval
f4 is step function with a discontinuity at 1
f5 is a hyperbola with vertical asymptote at 1
f6 has random values positive to left of 1, negative to right

of course these are not real problems. They just test how the
'good' solvers behave in bad circumstances where bisection is
really the best. A good solver should not be much worse than
bisection in such circumstance, while being faster for smooth
monotone sorts of functions.
"""

methods = [cc.bisect,cc.ridder,cc.brenth,cc.brentq]
mstrings = ['cc.bisect','cc.ridder','cc.brenth','cc.brentq']
functions = [f2,f3,f4,f5,f6]
fstrings = ['f2','f3','f4','f5','f6']

class TestBasic(TestCase) :
    def run_check(self, method, name):
        a = .5
        b = sqrt(3)
        for function, fname in zip(functions, fstrings):
            zero, r = method(function, a, b, xtol=0.1e-12, full_output=True)
            assert r.converged
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

    @dec.bench
    def test_run(self):
        a = .5
        b = sqrt(3)
        repeat = 2000

        print description

        print 'TESTING SPEED\n'
        print 'times in seconds for %d iterations \n'%repeat
        for i in range(len(functions)) :
            print 'function %s\n'%fstrings[i]
            func = functions[i]
            for j in range(len(methods)) :
                meth = methods[j]
                try:
                    t = measure("meth(func,a,b)",repeat)
                except:
                    print '%s : failed'%mstrings[j]
                else:
                    print '%s : %5.3f'%(mstrings[j],t)
            print '\n\n'

if __name__ == '__main__' :
    nose.run(argv=['', __file__])
