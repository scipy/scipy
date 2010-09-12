#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_array_almost_equal, assert_
from numpy import ones

from scipy.linalg import flapack, clapack


class TestFlapackSimple(TestCase):

    def test_gebal(self):
        a = [[1,2,3],[4,5,6],[7,8,9]]
        a1 = [[1,0,0,3e-4],
              [4,0,0,2e-3],
              [7,1,0,0],
              [0,1,0,0]]
        for p in 'sdzc':
            f = getattr(flapack,p+'gebal',None)
            if f is None: continue
            ba,lo,hi,pivscale,info = f(a)
            assert_(not info,`info`)
            assert_array_almost_equal(ba,a)
            assert_equal((lo,hi),(0,len(a[0])-1))
            assert_array_almost_equal(pivscale,ones(len(a)))

            ba,lo,hi,pivscale,info = f(a1,permute=1,scale=1)
            assert_(not info,`info`)
            #print a1
            #print ba,lo,hi,pivscale

    def test_gehrd(self):
        a = [[-149, -50,-154],
             [ 537, 180, 546],
             [ -27,  -9, -25]]
        for p in 'd':
            f = getattr(flapack,p+'gehrd',None)
            if f is None: continue
            ht,tau,info = f(a)
            assert_(not info,`info`)

class TestLapack(TestCase):

    def test_flapack(self):
        if hasattr(flapack,'empty_module'):
            #flapack module is empty
            pass

    def test_clapack(self):
        if hasattr(clapack,'empty_module'):
            #clapack module is empty
            pass


if __name__ == "__main__":
    run_module_suite()
