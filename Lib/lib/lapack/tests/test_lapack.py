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
from scipy_base import ones,dot,identity
set_package_path()
from lapack import flapack,clapack
restore_path()

class _test_ev:

    def check_syev(self,level=1,sym='sy',suffix=''):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        f = getattr(self.lapack,sym+'ev'+suffix)
        w,v,info=f(a)
        assert not info,`info`
        assert_array_almost_equal(w,exact_w)
        for i in range(3):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_syevd(self):
        self.check_syev(suffix='d')

    def check_heev(self):
        self.check_syev(sym='he')

    def check_heevd(self):
        self.check_syev(sym='he',suffix='d')

    def check_heev_complex(self,level=1,suffix=''):
        a= [[1,2-2j,3+7j],[2+2j,2,3],[3-7j,3,5]]
        exact_w=[-6.305141710654834,2.797880950890922,11.50726075976392]
        f = getattr(self.lapack,'heev'+suffix)
        w,v,info=f(a)
        assert not info,`info`
        assert_array_almost_equal(w,exact_w)
        for i in range(3):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_heevd_complex(self):
        self.check_heev_complex(suffix='d')

    def check_syevr(self,level=1,sym='sy'):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        f = getattr(self.lapack,sym+'evr')
        w,v,info = f(a)
        assert not info,`info`
        assert_array_almost_equal(w,exact_w)
        for i in range(3):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_heevr_complex(self,level=1):
        a= [[1,2-2j,3+7j],[2+2j,2,3],[3-7j,3,5]]
        exact_w=[-6.305141710654834,2.797880950890922,11.50726075976392]
        f = self.lapack.heevr
        w,v,info = f(a)
        assert not info,`info`
        assert_array_almost_equal(w,exact_w)
        for i in range(3):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_heevr(self):
        self.check_syevr(sym='he')

    def check_syevr_irange(self,level=1,sym='sy',irange=[0,2]):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        f = getattr(self.lapack,sym+'evr')
        w,v,info = f(a,irange=irange)
        assert not info,`info`
        rslice = slice(irange[0],irange[1]+1)
        m = irange[1] - irange[0] + 1
        assert_equal(len(w),m)
        assert_array_almost_equal(w,exact_w[rslice])
        for i in range(m):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_syevr_irange_low(self):
        self.check_syevr_irange(irange=[0,1])

    def check_syevr_irange_mid(self):
        self.check_syevr_irange(irange=[1,1])

    def check_syevr_irange_high(self):
        self.check_syevr_irange(irange=[1,2])

    def check_heevr_irange(self):
        self.check_syevr_irange(sym='he')

    def check_heevr_irange_low(self):
        self.check_syevr_irange(sym='he',irange=[0,1])

    def check_heevr_irange_high(self):
        self.check_syevr_irange(sym='he',irange=[1,2])

    def check_syevr_vrange(self,level=1,sym='sy',vrange=None):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        if vrange is None:
            vrange = [-1,10]
        ew = [value for value in exact_w if vrange[0]<value<=vrange[1]]
        f = getattr(self.lapack,sym+'evr')
        w,v,info = f(a,vrange=vrange)
        assert not info,`info`
        assert_array_almost_equal(w,ew)
        m = len(w)
        for i in range(m):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_syevr_vrange_low(self):
        self.check_syevr_vrange(vrange=[-1,1])

    def check_syevr_vrange_mid(self):
        self.check_syevr_vrange(vrange=[0,1])

    def check_syevr_vrange_high(self):
        self.check_syevr_vrange(vrange=[1,10])

    def check_heevr_vrange(self):
        self.check_syevr_vrange(sym='he')

    def check_heevr_vrange_low(self):
        self.check_syevr_vrange(sym='he',vrange=[-1,1])

    def check_heevr_vrange_high(self):
        self.check_syevr_vrange(sym='he',vrange=[1,10])

class _test_gev:

    def check_sygv(self,level=1,sym='sy',suffix='',itype=1):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        b = [[10,-1,1],[-1,8,-2],[1,-2,6]]
        f = getattr(self.lapack,sym+'gv'+suffix)
        w,v,info=f(a,b,itype=itype)
        assert not info,`info`
        for i in range(3):
            if itype==1:
                assert_array_almost_equal(dot(a,v[:,i]),w[i]*dot(b,v[:,i]),self.decimal)
            elif itype==2:
                assert_array_almost_equal(dot(a,dot(b,v[:,i])),w[i]*v[:,i],self.decimal)
            elif itype==3:
                assert_array_almost_equal(dot(b,dot(a,v[:,i])),w[i]*v[:,i],self.decimal-1)
            else:
                raise ValueError,`itype`

    def check_sygv_2(self): self.check_sygv(itype=2)

    def check_sygv_3(self): self.check_sygv(itype=3)

    def check_hegv(self): self.check_sygv(sym='he')

    def check_hegv_2(self): self.check_sygv(sym='he',itype=2)

    def check_hegv_3(self): self.check_sygv(sym='he',itype=3)

#class _test_ev: pass

class _test_lapack_simple(ScipyTestCase,_test_ev,_test_gev):

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

class PrefixWrapper:
    def __init__(self,module,prefix):
        self.module = module
        self.prefix = prefix
        self.__doc__ = module.__doc__

    def __getattr__(self, name):
        class A: pass
        a = getattr(self.module,self.prefix+name,getattr(self.module,name,A()))
        if isinstance(a,A):
            raise HideException,'%s has no attribute %r' % (self.module,name)
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
    class test_flapack_double(_test_lapack_simple):
        lapack = PrefixWrapper(flapack,'d')
        decimal = 12
    class test_flapack_float(_test_lapack_simple):
        lapack = PrefixWrapper(flapack,'s')
        decimal = 5
    class test_flapack_complex(_test_lapack_simple):
        lapack = PrefixWrapper(flapack,'c')
        decimal = 5
    class test_flapack_double_complex(_test_lapack_simple):
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
    class test_clapack_double(_test_lapack_simple):
        lapack = PrefixWrapper(clapack,'d')
        decimal = 12
    class test_clapack_float(_test_lapack_simple):
        lapack = PrefixWrapper(clapack,'s')
        decimal = 5
    class test_clapack_complex(_test_lapack_simple):
        lapack = PrefixWrapper(clapack,'c')
        decimal = 5
    class test_clapack_double_complex(_test_lapack_simple):
        lapack = PrefixWrapper(clapack,'z')
        decimal = 12

if __name__ == "__main__":
    ScipyTest().run()
