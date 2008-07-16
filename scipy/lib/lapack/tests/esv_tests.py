
from numpy.testing import *
from numpy import dot

class _test_ev(object):

    def check_syev(self,sym='sy',suffix=''):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        f = getattr(self.lapack,sym+'ev'+suffix)
        w,v,info=f(a)
        assert not info,`info`
        assert_array_almost_equal(w,exact_w)
        for i in range(3):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

    def check_syevd(self): self.check_syev(suffix='d')

    #def check_heev(self): self.check_syev(sym='he')

    #def check_heevd(self): self.check_syev(sym='he',suffix='d')

##    def check_heev_complex(self,suffix=''):
##        a= [[1,2-2j,3+7j],[2+2j,2,3],[3-7j,3,5]]
##        exact_w=[-6.305141710654834,2.797880950890922,11.50726075976392]
##        f = getattr(self.lapack,'heev'+suffix)
##        w,v,info=f(a)
##        assert not info,`info`
##        assert_array_almost_equal(w,exact_w)
##        for i in range(3):
##            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i],self.decimal)

    #def check_heevd_complex(self): self.check_heev_complex(suffix='d')

    def check_syevr(self,sym='sy'):
        a = [[1,2,3],[2,2,3],[3,3,6]]
        exact_w = [-0.6699243371851365,0.4876938861533345,9.182230451031804]
        f = getattr(self.lapack,sym+'evr')
        w,v,info = f(a)
        assert not info,`info`
        assert_array_almost_equal(w,exact_w)
        for i in range(3):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

##    def check_heevr_complex(self):
##        a= [[1,2-2j,3+7j],[2+2j,2,3],[3-7j,3,5]]
##        exact_w=[-6.305141710654834,2.797880950890922,11.50726075976392]
##        f = self.lapack.heevr
##        w,v,info = f(a)
##        assert not info,`info`
##        assert_array_almost_equal(w,exact_w)
##        for i in range(3):
##            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])

##    def check_heevr(self): self.check_syevr(sym='he')

    def check_syevr_irange(self,sym='sy',irange=[0,2]):
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

    def check_syevr_irange_low(self): self.check_syevr_irange(irange=[0,1])

    def check_syevr_irange_mid(self): self.check_syevr_irange(irange=[1,1])

    def check_syevr_irange_high(self): self.check_syevr_irange(irange=[1,2])

##    def check_heevr_irange(self): self.check_syevr_irange(sym='he')

##    def check_heevr_irange_low(self): self.check_syevr_irange(sym='he',irange=[0,1])

##    def check_heevr_irange_high(self): self.check_syevr_irange(sym='he',irange=[1,2])

    def check_syevr_vrange(self,sym='sy',vrange=None):
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

    def check_syevr_vrange_low(self): self.check_syevr_vrange(vrange=[-1,1])

    def check_syevr_vrange_mid(self): self.check_syevr_vrange(vrange=[0,1])

    def check_syevr_vrange_high(self): self.check_syevr_vrange(vrange=[1,10])

##    def check_heevr_vrange(self): self.check_syevr_vrange(sym='he')

##    def check_heevr_vrange_low(self): self.check_syevr_vrange(sym='he',vrange=[-1,1])

##    def check_heevr_vrange_high(self): self.check_syevr_vrange(sym='he',vrange=[1,10])
