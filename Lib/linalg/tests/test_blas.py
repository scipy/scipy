#!/usr/bin/env python
#
# Created by: Pearu Peterson, April 2002
#

__usage__ = """
Build linalg:
  python setup_linalg.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.linalg.test(<level>)'
Run tests if linalg is not installed:
  python tests/test_blas.py [<level>]
"""


from Numeric import arange, add, array
import math

import sys
from scipy_test.testing import *
set_package_path()
from linalg import fblas
from linalg import cblas
del sys.path[0]

class test_cblas1_simple(ScipyTestCase):

    def check_axpy(self):
        for p in 'sd':
            f = getattr(cblas,p+'axpy',None)
            if f is None: continue
            assert_array_almost_equal(f(5,[1,2,3],[2,-1,3]),[7,9,18])
        for p in 'cz':
            f = getattr(cblas,p+'axpy',None)
            if f is None: continue
            assert_array_almost_equal(f(5,[1,2j,3],[2,-1,3]),[7,10j-1,18])

class test_fblas1_simple(ScipyTestCase):

    def check_axpy(self):
        for p in 'sd':
            f = getattr(fblas,p+'axpy',None)
            if f is None: continue
            assert_array_almost_equal(f([1,2,3],[2,-1,3],a=5),[7,9,18])
        for p in 'cz':
            f = getattr(fblas,p+'axpy',None)
            if f is None: continue
            assert_array_almost_equal(f([1,2j,3],[2,-1,3],a=5),[7,10j-1,18])
    def check_copy(self):
        for p in 'sd':
            f = getattr(fblas,p+'copy',None)
            if f is None: continue
            assert_array_almost_equal(f([3,4,5],[8]*3),[3,4,5])
        for p in 'cz':
            f = getattr(fblas,p+'copy',None)
            if f is None: continue
            assert_array_almost_equal(f([3,4j,5+3j],[8]*3),[3,4j,5+3j])
    def check_asum(self):
        for p in 'sd':
            f = getattr(fblas,p+'asum',None)
            if f is None: continue
            assert_almost_equal(f([3,-4,5]),12)
        for p in ['sc','dz']:
            f = getattr(fblas,p+'asum',None)
            if f is None: continue
            assert_almost_equal(f([3j,-4,3-4j]),14)
    def check_dot(self):
        for p in 'sd':
            f = getattr(fblas,p+'dot',None)
            if f is None: continue
            assert_almost_equal(f([3,-4,5],[2,5,1]),-9)
        for p in 'cz':
            f = getattr(fblas,p+'dotu',None)
            if f is None: continue
            assert_almost_equal(f([3j,-4,3-4j],[2,3,1]),-9+2j)
            f = getattr(fblas,p+'dotc')
            assert_almost_equal(f([3j,-4,3-4j],[2,3j,1]),3-14j)
    def check_nrm2(self):
        for p in 'sd':
            f = getattr(fblas,p+'nrm2',None)
            if f is None: continue
            assert_almost_equal(f([3,-4,5]),math.sqrt(50))
        for p in ['sc','dz']:
            f = getattr(fblas,p+'nrm2',None)
            if f is None: continue
            assert_almost_equal(f([3j,-4,3-4j]),math.sqrt(50))
    def check_scal(self):
        for p in 'sd':
            f = getattr(fblas,p+'scal',None)
            if f is None: continue
            assert_array_almost_equal(f(2,[3,-4,5]),[6,-8,10])
        for p in 'cz':
            f = getattr(fblas,p+'scal',None)
            if f is None: continue
            assert_array_almost_equal(f(3j,[3j,-4,3-4j]),[-9,-12j,12+9j])
        for p in ['cs','zd']:
            f = getattr(fblas,p+'scal',None)
            if f is None: continue
            assert_array_almost_equal(f(3,[3j,-4,3-4j]),[9j,-12,9-12j])
    def check_swap(self):
        for p in 'sd':
            f = getattr(fblas,p+'swap',None)
            if f is None: continue
            x,y = [2,3,1],[-2,3,7]
            x1,y1 = f(x,y)
            assert_array_almost_equal(x1,y)
            assert_array_almost_equal(y1,x)
        for p in 'cz':
            f = getattr(fblas,p+'swap',None)
            if f is None: continue
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

class test_fblas2_simple(ScipyTestCase):

    def check_gemv(self):
        for p in 'sd':
            f = getattr(fblas,p+'gemv',None)
            if f is None: continue
            assert_array_almost_equal(f(3,[[3]],[-4]),[-36])
            assert_array_almost_equal(f(3,[[3]],[-4],3,[5]),[-21])
        for p in 'cz':
            f = getattr(fblas,p+'gemv',None)
            if f is None: continue
            assert_array_almost_equal(f(3j,[[3-4j]],[-4]),[-48-36j])
            assert_array_almost_equal(f(3j,[[3-4j]],[-4],3,[5j]),[-48-21j])

    def check_ger(self):

        for p in 'sd':
            f = getattr(fblas,p+'ger',None)
            if f is None: continue
            assert_array_almost_equal(f(1,[1,
                                           2],[3,4]),[[3,4],[6,8]])
            assert_array_almost_equal(f(2,[1,
                                           2,
                                           3],[3,4]),[[6,8],[12,16],[18,24]])

            assert_array_almost_equal(f(1,[1,
                                           2],[3,4],
                                        a=[[1,2],[3,4]]
                                        ),[[4,6],[9,12]])

        for p in 'cz':
            f = getattr(fblas,p+'geru',None)
            if f is None: continue
            assert_array_almost_equal(f(1,[1j,
                                           2],[3,4]),[[3j,4j],[6,8]])
            assert_array_almost_equal(f(-2,[1j,
                                           2j,
                                           3j],[3j,4j]),[[6,8],[12,16],[18,24]])

        for p in 'cz':
            f = getattr(fblas,p+'gerc',None)
            if f is None: continue
            assert_array_almost_equal(f(1,[1j,
                                           2],[3,4]),[[3j,4j],[6,8]])
            assert_array_almost_equal(f(2,[1j,
                                           2j,
                                           3j],[3j,4j]),[[6,8],[12,16],[18,24]])

class test_fblas3_simple(ScipyTestCase):

    def check_gemm(self):
        for p in 'sd':
            f = getattr(fblas,p+'gemm',None)
            if f is None: continue
            assert_array_almost_equal(f(3,[3],[-4]),[-36])
            assert_array_almost_equal(f(3,[3],[-4],3,[5]),[-21])
        for p in 'cz':
            f = getattr(fblas,p+'gemm',None)
            if f is None: continue
            assert_array_almost_equal(f(3j,[3-4j],[-4]),[-48-36j])
            assert_array_almost_equal(f(3j,[3-4j],[-4],3,[5j]),[-48-21j])

class test_blas(ScipyTestCase):

    def check_fblas(self):
        if hasattr(fblas,'empty_module'):
            print """
****************************************************************
WARNING: fblas module is empty.
-----------
See scipy/INSTALL.txt for troubleshooting.
****************************************************************
"""
    def check_cblas(self):
        if hasattr(cblas,'empty_module'):
            print """
****************************************************************
WARNING: cblas module is empty
-----------
See scipy/INSTALL.txt for troubleshooting.
Notes:
* If atlas library is not found by scipy/system_info.py,
  then scipy uses fblas instead of cblas.
****************************************************************
"""

if __name__ == "__main__":
    ScipyTest('linalg.blas').run()
