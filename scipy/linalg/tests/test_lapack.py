#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_array_almost_equal, assert_

import numpy as np

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
            assert_array_almost_equal(pivscale, np.ones(len(a)))

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
            
    def test_trsyl(self):
        a = [[1, 2], [3, 4]]
        b = [[5, 6], [7, 8]]
        c = [[9, 10], [11, 12]]
        
        x_expected = [[0.5, 0.66667],
                      [0.66667, 0.5]]
                      
        x_expected_t = [[0.5, 0.58333],
                        [0.83333, 0.41667]]
                        
        x_expected_m = [[-1.125, -1.0],
                        [-1.25, -1.875]]
        
        # Test single and double implementations, including most
        # of the options
        for p in 'ds':
            f = getattr(flapack, p+'trsyl', None)
                
            x, scale, info = f(a, b, c)
            assert_array_almost_equal(x, x_expected, decimal=4)
            assert_equal(scale,1.0)
            
            x, scale, info = f(a, b, c, trana='T', tranb='T')
            assert_array_almost_equal(x, x_expected_t, decimal=4)
            
            x, scale, info = f(a, b, c, isgn=-1)
            assert_array_almost_equal(x, x_expected_m, decimal=4)
        
        # Test complex operation (a bit simpler...)
        ac = [[0+2.j,],]
        bc = [[0+4.j,],]
        cc = [[0+12.j,],]
        x_expected_c = [[2.0+0.j,],]
        
        for p in 'zc':
            f = getattr(flapack, p+'trsyl', None)
            
            x, scale, info = f(ac, bc, cc)
            assert_array_almost_equal(x, x_expected_c, decimal=4)
            
            x, scale, info = f(ac, bc, cc, trana='C', tranb='C')
            assert_array_almost_equal(x, -1.0*np.asarray(x_expected_c), decimal=4)

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
