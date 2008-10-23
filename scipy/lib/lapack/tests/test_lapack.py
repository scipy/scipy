#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#
'''
This file adapted for nose tests 1/1/08

Note that the conversion is not very complete.

This and the included files deliberately use "check_" as the test
method names.  There are no subclasses of TestCase.  Thus nose will
pick up nothing but the final test_all_lapack generator function.
This does the work of collecting the test methods and checking if they
can be run (see the isrunnable method).
'''

from numpy.testing import *
from numpy import ones

from scipy.lib.lapack import flapack, clapack

#sys.path.insert(0, os.path.split(__file__))
from gesv_tests import _test_gev
from esv_tests import _test_ev
#del sys.path[0]

#class _test_ev: pass

class _TestLapack( _test_ev,
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
        ''' Return True if required routines for check method present in module '''
        l = mthname.split('_')
        if len(l)>1 and l[0]=='check':
            return hasattr(self.lapack,l[1])
        return 2

class PrefixWrapper(object):
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
    class TestFlapackDouble(_TestLapack):
        lapack = PrefixWrapper(flapack,'d')
        decimal = 12
    class TestFlapackFloat(_TestLapack):
        lapack = PrefixWrapper(flapack,'s')
        decimal = 5
    class TestFlapackComplex(_TestLapack):
        lapack = PrefixWrapper(flapack,'c')
        decimal = 5
    class TestFlapackDoubleComplex(_TestLapack):
        lapack = PrefixWrapper(flapack,'z')
        decimal = 12

if hasattr(clapack,'empty_module') or clapack is flapack:
    print """
****************************************************************
WARNING: clapack module is empty
-----------
See scipy/INSTALL.txt for troubleshooting.
Notes:
* If atlas library is not found by numpy/distutils/system_info.py,
  then scipy uses flapack instead of clapack.
****************************************************************
"""
else:
    class TestClapackDouble(_TestLapack):
        lapack = PrefixWrapper(clapack,'d')
        decimal = 12
    class TestClapackFloat(_TestLapack):
        lapack = PrefixWrapper(clapack,'s')
        decimal = 5
    class TestClapackComplex(_TestLapack):
        lapack = PrefixWrapper(clapack,'c')
        decimal = 5
    class TestClapackDoubleComplex(_TestLapack):
        lapack = PrefixWrapper(clapack,'z')
        decimal = 12

# Collect test classes and methods with generator
# This is a moderate hack replicating some obscure numpy testing
# functionality for use with nose

def test_all_lapack():
    methods = []
    for name, value in globals().items():
        if not (name.startswith('Test')
                and issubclass(value, _TestLapack)):
            continue
        o = value()
        methods += [getattr(o, n) for n in dir(o) if o.isrunnable(n) is True]
    for method in methods:
        yield (method, )
