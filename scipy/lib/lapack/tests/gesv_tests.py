
from numpy.testing import *
from numpy import dot

class _test_gev(object):

    def check_sygv(self,sym='sy',suffix='',itype=1):
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

##    def check_hegv(self): self.check_sygv(sym='he')

##    def check_hegv_2(self): self.check_sygv(sym='he',itype=2)

##    def check_hegv_3(self): self.check_sygv(sym='he',itype=3)

    def check_sygvd(self): self.check_sygv(suffix='d')

    def check_sygvd_2(self): self.check_sygv(suffix='d',itype=2)

    def check_sygvd_3(self): self.check_sygv(suffix='d',itype=3)

##    def check_hegvd(self): self.check_sygv(sym='he',suffix='d')

##    def check_hegvd_2(self): self.check_sygv(sym='he',suffix='d',itype=2)

##    def check_hegvd_3(self): self.check_sygv(sym='he',suffix='d',itype=3)
