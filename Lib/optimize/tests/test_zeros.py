#!/usr/local/lib/python

import sys,time
import scipy.optimize.zeros as cc
from math import sin,sqrt,log
from whrandom import random

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

description =\
"""
f2 is a symmetric parabola, x**2 - 1
f3 is a quartic polynomial with large hump in interval
f4 is step function with a discontinuity at 1
f5 is a hyperbola with vertical asymptote at 1
f6 has random values positive to left of 1 , negative to right

of course these are not real problems. They just test how the
'good' solvers behave in bad circumstances where bisection is
really the best. A good solver should not be much worse than
bisection in such circumstance, while being faster for smooth
monotone sorts of functions.
"""

methods = [cc.bisect,cc.ridder,cc.brenth,cc.brent]
mstrings = ['cc.bisect','cc.ridder','cc.brenth','cc.brent']
functions = [f2,f3,f4,f5,f6]
fstrings = ['f2','f3','f4','f5','f6']

from scipy_base.testing import ScipyTestCase
import unittest

class test_basic(ScipyTestCase) :

    def check_run(self) :
        a = .5
        b = sqrt(3)
        repeat = 2000

        
        print 'TESTING CONVERGENCE\n'
        print 'zero should be 1\n'
        for i in range(len(functions)) :
            print 'function %s\n'%fstrings[i]
            for j in range(len(methods)) :
                try :
                    zero = methods[j](functions[i],a,b)
                except :
                    print '%s : failed'%mstrings[j]
                else:
                    print '%s : %21.16f'%(mstrings[j],zero)
            print '\n\n'

    def bench_run(self):
        a = .5
        b = sqrt(3)
        repeat = 2000

        print '%s\n',description

        print 'TESTING SPEED\n'
        print 'times in seconds for %d iterations \n'%repeat
        for i in range(len(functions)) :
            print 'function %s\n'%fstrings[i]
            func = functions[i]
            for j in range(len(methods)) :
                meth = methods[j]
                try :
                    time = self.measure("meth(func,a,b)",repeat)
                except :
                    print '%s : failed'%mstrings[j]
                else:
                    print '%s : %5.3f'%(mstrings[j],time)
            print '\n\n'


def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_basic,'check_'))
    if level > 5:
        suites.append( unittest.makeSuite(test_basic,'bench_'))

    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level=level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner

if __name__ == '__main__' :
    if len(sys.argv)>1:
        level = eval(sys.argv[1])
    else:
        level = 10    
    test(level)
