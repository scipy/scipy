#!/usr/bin/env python
# Usage:
#   In the parent directory run
#     python setup_linalg.py build --build-platlib=.
#     python -c 'import linalg;linalg.test(1)'
#
# Created by: Pearu Peterson, April 2002
#

from Numeric import arange, add, array
import math

from scipy_base.testing import assert_array_almost_equal, assert_equal
from scipy_base.testing import assert_almost_equal, assert_array_equal
from scipy_base.testing import ScipyTestCase
import unittest


import os,sys
d = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,d)
import cblas
import fblas
del sys.path[0]


class test_blas1_simple(ScipyTestCase):

    def check_axpy(self):
        for p in 'sd':
            f = getattr(cblas,p+'axpy')
            assert_array_almost_equal(f(5,[1,2,3],[2,-1,3]),[7,9,18])
            f = getattr(fblas,p+'axpy')
            assert_array_almost_equal(f(5,[1,2,3],[2,-1,3]),[7,9,18])
        for p in 'cz':
            f = getattr(cblas,p+'axpy')
            assert_array_almost_equal(f(5,[1,2j,3],[2,-1,3]),[7,10j-1,18])
            f = getattr(fblas,p+'axpy')
            assert_array_almost_equal(f(5,[1,2j,3],[2,-1,3]),[7,10j-1,18])
    def check_copy(self):
        for p in 'sd':
            f = getattr(fblas,p+'copy')
            assert_array_almost_equal(f([3,4,5],[8]*3),[3,4,5])
        for p in 'cz':
            f = getattr(fblas,p+'copy')
            assert_array_almost_equal(f([3,4j,5+3j],[8]*3),[3,4j,5+3j])
    def check_asum(self):
        for p in 'sd':
            f = getattr(fblas,p+'asum')
            assert_almost_equal(f([3,-4,5]),12)
        for p in ['sc','dz']:
            f = getattr(fblas,p+'asum')
            assert_almost_equal(f([3j,-4,3-4j]),14)
    def check_dot(self):
        for p in 'sd':
            f = getattr(fblas,p+'dot')
            assert_almost_equal(f([3,-4,5],[2,5,1]),-9)
        for p in 'cz':
            f = getattr(fblas,p+'dotu')
            assert_almost_equal(f([3j,-4,3-4j],[2,3,1]),-9+2j)
            f = getattr(fblas,p+'dotc')
            assert_almost_equal(f([3j,-4,3-4j],[2,3j,1]),3-14j)
    def check_nrm2(self):
        for p in 'sd':
            f = getattr(fblas,p+'nrm2')
            assert_almost_equal(f([3,-4,5]),math.sqrt(50))
        for p in ['sc','dz']:
            f = getattr(fblas,p+'nrm2')
            assert_almost_equal(f([3j,-4,3-4j]),math.sqrt(50))
    def check_scal(self):
        for p in 'sd':
            f = getattr(fblas,p+'scal')
            assert_array_almost_equal(f(2,[3,-4,5]),[6,-8,10])
        for p in 'cz':
            f = getattr(fblas,p+'scal')
            assert_array_almost_equal(f(3j,[3j,-4,3-4j]),[-9,-12j,12+9j])
        for p in ['cs','zd']:
            f = getattr(fblas,p+'scal')
            assert_array_almost_equal(f(3,[3j,-4,3-4j]),[9j,-12,9-12j])
    def check_swap(self):
        for p in 'sd':
            f = getattr(fblas,p+'swap')
            x,y = [2,3,1],[-2,3,7]
            x1,y1 = f(x,y)
            assert_array_almost_equal(x1,y)
            assert_array_almost_equal(y1,x)
        for p in 'cz':
            f = getattr(fblas,p+'swap')
            x,y = [2,3j,1],[-2,3,7-3j]
            x1,y1 = f(x,y)
            assert_array_almost_equal(x1,y)
            assert_array_almost_equal(y1,x)
    def check_amax(self):
        for p in 'sd':
            f = getattr(fblas,'i'+p+'amax')
            assert_equal(f([-2,4,3]),1)
        for p in 'cz':
            f = getattr(fblas,'i'+p+'amax')
            assert_equal(f([-5,4+3j,6]),1)
    #XXX: need tests for rot,rotm,rotg,rotmg

class test_blas2_simple(ScipyTestCase):

    def check_gemv(self):
        for p in 'sd':
            f = getattr(fblas,p+'gemv')
            assert_array_almost_equal(f(3,[[3]],[-4]),[-36])
            assert_array_almost_equal(f(3,[[3]],[-4],3,[5]),[-21])
        for p in 'cz':
            f = getattr(fblas,p+'gemv')
            assert_array_almost_equal(f(3j,[[3-4j]],[-4]),[-48-36j])
            assert_array_almost_equal(f(3j,[[3-4j]],[-4],3,[5j]),[-48-21j])


class test_blas3_simple(ScipyTestCase):

    def check_gemm(self):
        for p in 'sd':
            f = getattr(fblas,p+'gemm')
            assert_array_almost_equal(f(3,[3],[-4]),[-36])
            assert_array_almost_equal(f(3,[3],[-4],3,[5]),[-21])
        for p in 'cz':
            f = getattr(fblas,p+'gemm')
            assert_array_almost_equal(f(3j,[3-4j],[-4]),[-48-36j])
            assert_array_almost_equal(f(3j,[3-4j],[-4],3,[5j]),[-48-21j])

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_blas1_simple,'check_') )
        suites.append( unittest.makeSuite(test_blas2_simple,'check_') )
        suites.append( unittest.makeSuite(test_blas3_simple,'check_') )
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner
