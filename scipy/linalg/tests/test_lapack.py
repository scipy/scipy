#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_array_almost_equal, assert_

import numpy as np

from scipy.linalg import flapack, clapack
from scipy.linalg.lapack import get_lapack_funcs


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
        a = np.array([[1, 2], [0, 4]])
        b = np.array([[5, 6], [0, 8]])
        c = np.array([[9, 10], [11, 12]])
        trans = 'T'

        # Test single and double implementations, including most
        # of the options
        for dtype in 'fdFD':
            a1, b1, c1 = a.astype(dtype), b.astype(dtype), c.astype(dtype)
            trsyl, = get_lapack_funcs(('trsyl',), (a1,))
            if dtype.isupper(): # is complex dtype
                a1[0] += 1j
                trans = 'C'

            x, scale, info = trsyl(a1, b1, c1)
            assert_array_almost_equal(np.dot(a1, x) + np.dot(x, b1), scale * c1)

            x, scale, info = trsyl(a1, b1, c1, trana=trans, tranb=trans)
            assert_array_almost_equal(np.dot(a1.conjugate().T, x) + np.dot(x, b1.conjugate().T),
                scale * c1, decimal=4)

            x, scale, info = trsyl(a1, b1, c1, isgn=-1)
            assert_array_almost_equal(np.dot(a1, x) - np.dot(x, b1), scale * c1, decimal=4)

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
