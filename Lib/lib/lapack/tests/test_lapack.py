#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#

__usage__ = """
Build lapack:
  python setup_lapack.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.lib.lapack.test(<level>)'
Run tests if lapack is not installed:
  python tests/test_lapack.py [<level>]
"""

import sys
from scipy_test.testing import *
from scipy_base import ones,dot
set_package_path()
from lapack import flapack,clapack
restore_path()

class _test_lapack_simple(ScipyTestCase):

    def check_syev(self,level=1,suffix=''):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        for p in 'sd':
            f = getattr(self.lapack,p+'syev'+suffix,None)
            if f is None: continue
            w,v,info=f(a)
            assert not info,`info`
            assert_array_almost_equal(w,exact_w)
            for i in range(3):
                assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_syevd(self):
        self.check_syev(suffix='d')

    def check_heev(self,level=1,suffix=''):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        for p in 'cz':
            f = getattr(self.lapack,p+'heev'+suffix,None)
            if f is None: continue
            w,v,info=f(a)
            assert not info,`info`
            assert_array_almost_equal(w,exact_w)
            for i in range(3):
                assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_heevd(self):
        self.check_heev(suffix='d')

    def check_gebal(self):
        a = [[1,2,3],[4,5,6],[7,8,9]]
        a1 = [[1,0,0,3e-4],
              [4,0,0,2e-3],
              [7,1,0,0],
              [0,1,0,0]]
        for p in 'sdzc':
            f = getattr(self.lapack,p+'gebal',None)
            if f is None: continue
            ba,lo,hi,pivscale,info = f(a)
            assert not info,`info`
            assert_array_almost_equal(ba,a)
            assert_equal((lo,hi),(0,len(a[0])-1))
            assert_array_almost_equal(pivscale,ones(len(a)))

            ba,lo,hi,pivscale,info = f(a1,permute=1,scale=1)
            assert not info,`info`

    def check_gehrd(self):
        a = [[-149, -50,-154],
             [ 537, 180, 546],
             [ -27,  -9, -25]]
        for p in 'sdzc':
            f = getattr(self.lapack,p+'gehrd',None)
            if f is None: continue
            ht,tau,info = f(a)
            assert not info,`info`

if hasattr(flapack,'empty_module'):
    print """
****************************************************************
WARNING: flapack module is empty
-----------
See scipy/INSTALL.txt for troubleshooting.
****************************************************************
"""
else:
    class test_flapack_simple(_test_lapack_simple):
        lapack = flapack

if hasattr(clapack,'empty_module'):
    print """
****************************************************************
WARNING: clapack module is empty
-----------
See scipy/INSTALL.txt for troubleshooting.
Notes:
* If atlas library is not found by scipy/system_info.py,
  then scipy uses flapack instead of clapack.
****************************************************************
"""
    class test_clapack_simple(_test_lapack_simple):
        lapack = clapack


if __name__ == "__main__":
    ScipyTest().run()
