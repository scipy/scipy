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
from scipy.base import *

set_package_path()
from lapack import flapack,clapack
restore_path()

set_local_path()
from gesv_tests import _test_gev
from esv_tests import _test_ev
restore_path()

#class _test_ev: pass

class _test_lapack(ScipyTestCase,
                   _test_ev,
                   _test_gev):

    def check_gebal(self):
        a = [[1,2,3],[4,5,6],[7,8,9]]
        a1 = [[1,0,0,3e-4],
              [4,0,0,2e-3],
              [7,1,0,0],
              [0,1,0,0]]
        f = self.lapack.gebal

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
        f = self.lapack.gehrd
        ht,tau,info = f(a)
        assert not info,`info`

    def isrunnable(self,mthname):
        l = mthname.split('_')
        if len(l)>1 and l[0]=='check':
            return hasattr(self.lapack,l[1])
        return 2

class PrefixWrapper:
    def __init__(self,module,prefix):
        self.module = module
        self.prefix = prefix
        self.__doc__ = module.__doc__

    def __getattr__(self, name):
        class A: pass
        a = getattr(self.module,self.prefix+name,getattr(self.module,name,A()))
        if isinstance(a,A):
            raise AttributeError,'%s has no attribute %r' % (self.module,name)
        return a

if hasattr(flapack,'empty_module'):
    print """
****************************************************************
WARNING: flapack module is empty
-----------
See scipy/INSTALL.txt for troubleshooting.
****************************************************************
"""
else:
    class test_flapack_double(_test_lapack):
        lapack = PrefixWrapper(flapack,'d')
        decimal = 12
    class test_flapack_float(_test_lapack):
        lapack = PrefixWrapper(flapack,'s')
        decimal = 5
    class test_flapack_complex(_test_lapack):
        lapack = PrefixWrapper(flapack,'c')
        decimal = 5
    class test_flapack_double_complex(_test_lapack):
        lapack = PrefixWrapper(flapack,'z')
        decimal = 12

if hasattr(clapack,'empty_module') or clapack is flapack:
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
else:
    class test_clapack_double(_test_lapack):
        lapack = PrefixWrapper(clapack,'d')
        decimal = 12
    class test_clapack_float(_test_lapack):
        lapack = PrefixWrapper(clapack,'s')
        decimal = 5
    class test_clapack_complex(_test_lapack):
        lapack = PrefixWrapper(clapack,'c')
        decimal = 5
    class test_clapack_double_complex(_test_lapack):
        lapack = PrefixWrapper(clapack,'z')
        decimal = 12

if __name__ == "__main__":
    ScipyTest().run()
