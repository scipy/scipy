#this program corresponds to special.py

### Means test is not done yet
#E   Means test is giving error (E)
#F   Means test is failing (F)
#EF  Means test is giving error and Failing
#!   Means test is segfaulting

###  test_besselpoly
###  test_jnjnp_zeros
###  test_mathieu_a
###  test_mathieu_even_coef
###  test_mathieu_odd_coef
###  test_modfresnelp
###  test_modfresnelm
#    test_pbdv_seq
###  test_pbvv_seq
###  test_sph_harm
#    test_sph_in
#    test_sph_jn
#    test_sph_kn

#With Py2.3:
#!   test_sh_chebyt
#!   test_pbdv_seq
#F   test_sph_kn
#F   test_sph_jn
#F   test_sph_in
#F   test_sph_jn

import os
import sys
#import fpformat
import unittest
import scipy_base.limits as limits
from scipy_base import *

from scipy_test.testing import *
set_package_path()
from special import *
del sys.path[0]


class test_cephes(ScipyTestCase):
    def check_airy(self):
        cephes.airy(0)
    def check_airye(self):
        cephes.airye(0)
    def check_bdtr(self):
        assert_equal(cephes.bdtr(1,1,0.5),1.0)
    def check_bdtri(self):
        assert_equal(cephes.bdtri(1,3,0.5),0.5)
    def check_bdtrc(self):
        assert_equal(cephes.bdtrc(1,3,0.5),0.5)
    def check_bdtrin(self):
        assert_equal(cephes.bdtrin(1,0,1),5.0)
    def check_bdtrik(self):
        cephes.bdtrik(1,3,0.5)

    def check_bei(self):
        assert_equal(cephes.bei(0),0.0)
    def check_beip(self):
        assert_equal(cephes.beip(0),0.0)
    def check_ber(self):
        assert_equal(cephes.ber(0),1.0)
    def check_berp(self):
        assert_equal(cephes.berp(0),0.0)

    def check_besselpoly(self):
        assert_equal(cephes.besselpoly(0,0,0),1.0)

    def check_beta(self):
        assert_equal(cephes.beta(1,1),1.0)
    def check_betainc(self):
        assert_equal(cephes.betainc(1,1,1),1.0)
    def check_betaln(self):
        assert_equal(cephes.betaln(1,1),0.0)
    def check_betaincinv(self):
        assert_equal(cephes.betaincinv(1,1,1),1.0)

    def check_btdtr(self):
        assert_equal(cephes.btdtr(1,1,1),1.0)
    def check_btdtri(self):
        assert_equal(cephes.btdtri(1,1,1),1.0)
    def check_btdtria(self):
        assert_equal(cephes.btdtria(1,1,1),5.0)
    def check_btdtrib(self):
        assert_equal(cephes.btdtrib(1,1,1),5.0)

    def check_cbrt(self):
        assert_approx_equal(cephes.cbrt(1),1.0)

    def check_chdtr(self):
        assert_equal(cephes.chdtr(1,0),0.0)
    def check_chdtrc(self):
        assert_equal(cephes.chdtrc(1,0),1.0)
    def check_chdtri(self):
        assert_equal(cephes.chdtri(1,1),0.0)
    def check_chdtriv(self):
        assert_equal(cephes.chdtriv(0,0),5.0)

    def check_chndtr(self):
        assert_equal(cephes.chndtr(0,1,0),0.0)
    def check_chndtridf(self):
        assert_equal(cephes.chndtridf(0,0,1),5.0)
    def check_chndtrinc(self):
        assert_equal(cephes.chndtrinc(0,1,0),5.0)
    def check_chndtrix(self):
        assert_equal(cephes.chndtrix(0,1,0),0.0)

    def check_cosdg(self):
        assert_equal(cephes.cosdg(0),1.0)
    def check_cosm1(self):
        assert_equal(cephes.cosm1(0),0.0)
    def check_cotdg(self):
        assert_almost_equal(cephes.cotdg(45),1.0)

    def check_dawsn(self):
        assert_equal(cephes.dawsn(0),0.0)

    def check_ellipe(self):
        assert_equal(cephes.ellipe(1),1.0)
    def check_ellipeinc(self):
        assert_equal(cephes.ellipeinc(0,1),0.0)
    def check_ellipj(self):
        cephes.ellipj(0,1)
    def check_ellipk(self):
        cephes.ellipk(0)#==pi/2
    def check_ellipkinc(self):
        assert_equal(cephes.ellipkinc(0,0),0.0)

    def check_erf(self):
        assert_equal(cephes.erf(0),0.0)
    def check_erfc(self):
        assert_equal(cephes.erfc(0),1.0)

    def check_exp1(self):
        cephes.exp1(1)
    def check_expi(self):
        cephes.expi(1)
    def check_expn(self):
        cephes.expn(1,1)

    def check_exp10(self):
        assert_approx_equal(cephes.exp10(2),100.0)
    def check_exp2(self):
        assert_equal(cephes.exp2(2),4.0)
    def check_expm1(self):
        assert_equal(cephes.expm1(0),0.0)

    def check_fdtr(self):
        assert_equal(cephes.fdtr(1,1,0),0.0)
    def check_fdtrc(self):
        assert_equal(cephes.fdtrc(1,1,0),1.0)
    def check_fdtri(self):
        cephes.fdtri(1,1,0.5)
    def check_fdtridfd(self):
        assert_equal(cephes.fdtridfd(1,0,0),5.0)

    def check_fresnel(self):
        assert_equal(cephes.fresnel(0),(0.0,0.0))
        
    def check_gamma(self):
        assert_equal(cephes.gamma(5),24.0)
    def check_gammainc(self):
        assert_equal(cephes.gammainc(5,0),0.0)
    def check_gammaincc(self):
        assert_equal(cephes.gammaincc(5,0),1.0)
    def check_gammainccinv(self):
        assert_equal(cephes.gammainccinv(5,1),0.0)
    def check_gammaln(self):
        cephes.gammaln(10)

    def check_gdtr(self):
        assert_equal(cephes.gdtr(1,1,0),0.0)
    def check_gdtrc(self):
        assert_equal(cephes.gdtrc(1,1,0),1.0)
    def check_gdtria(self):
        assert_equal(cephes.gdtria(0,1,1),0.0)
    def check_gdtrib(self):
        cephes.gdtrib(1,0,1)
        #assert_equal(cephes.gdtrib(1,0,1),5.0)
    def check_gdtrix(self):
        cephes.gdtrix(1,1,.1)

    def check_hankel1(self):
        cephes.hankel1(1,1)
    def check_hankel1e(self):
        cephes.hankel1e(1,1)
    def check_hankel2(self):
        cephes.hankel2(1,1)
    def check_hankel2e(self):
        cephes.hankel2e(1,1)

    def check_hyp1f1(self):
        cephes.hyp1f1(1,1,1)
    def check_hyp1f2(self):
        cephes.hyp1f2(1,1,1,1)
    def check_hyp2f0(self):
        cephes.hyp2f0(1,1,1,1)
    def check_hyp2f1(self):
        assert_equal(cephes.hyp2f1(1,1,1,0),1.0)
    def check_hyp3f0(self):
        assert_equal(cephes.hyp3f0(1,1,1,0),(1.0,0.0))
    def check_hyperu(self):
        assert_equal(cephes.hyperu(0,1,1),1.0)

    def check_i0(self):
        assert_equal(cephes.i0(0),1.0)
    def check_i0e(self):
        assert_equal(cephes.i0e(0),1.0)
    def check_i1(self):
        assert_equal(cephes.i1(0),0.0)
    def check_i1e(self):
        assert_equal(cephes.i1e(0),0.0)

    def check_it2i0k0(self):
        cephes.it2i0k0(1)
    def check_it2j0y0(self):
        cephes.it2j0y0(1)
    def check_it2struve0(self):
        cephes.it2struve0(1)
    def check_itairy(self):
        cephes.itairy(1)
    def check_iti0k0(self):
        assert_equal(cephes.iti0k0(0),(0.0,0.0))
    def check_itj0y0(self):
        assert_equal(cephes.itj0y0(0),(0.0,0.0))
    def check_itmodstruve0(self):
        assert_equal(cephes.itmodstruve0(0),0.0)
    def check_itstruve0(self):
        assert_equal(cephes.itstruve0(0),0.0)
    def check_iv(self):
        assert_equal(cephes.iv(1,0),0.0)
    def _check_ive(self):
        assert_equal(cephes.ive(1,0),0.0)

    def check_j0(self):
        assert_equal(cephes.j0(0),1.0)
    def check_j1(self):
        assert_equal(cephes.j1(0),0.0)
    def check_jn(self):
        assert_equal(cephes.jn(0,0),1.0)
    def check_jv(self):
        assert_equal(cephes.jv(0,0),1.0)
    def _check_jve(self):
        assert_equal(cephes.jve(0,0),1.0)

    def check_k0(self):
        cephes.k0(2)
    def check_k0e(self):
        cephes.k0e(2)
    def check_k1(self):
        cephes.k1(2)
    def check_k1e(self):
        cephes.k1e(2)
    def check_kei(self):
        cephes.kei(2)
    def check_keip(self):
        assert_equal(cephes.keip(0),0.0)
    def check_ker(self):
        cephes.ker(2)
    def check_kerp(self):
        cephes.kerp(2)
    def _check_kelvin(self):
        cephes.kelvin(2)
    def check_kn(self):
        cephes.kn(1,1)

    def check_kolmogi(self):
        assert_equal(cephes.kolmogi(1),0.0)
    def check_kolmogorov(self):
        assert_equal(cephes.kolmogorov(0),1.0)

    def _check_kv(self):
        cephes.kv(1,1)
    def _check_kve(self):
        cephes.kve(1,1)
    def check_log1p(self):
        assert_equal(cephes.log1p(0),0.0)
    def check_lpmv(self):
        assert_equal(cephes.lpmv(0,0,1),1.0)

    def check_mathieu_a(self):
        assert_equal(cephes.mathieu_a(1,0),1.0)
    def check_mathieu_b(self):
        assert_equal(cephes.mathieu_b(1,0),1.0)
    def check_mathieu_cem(self):
        assert_equal(cephes.mathieu_cem(1,0,0),(1.0,0.0))
    def check_mathieu_modcem1(self):
        assert_equal(cephes.mathieu_modcem1(1,0,0),(0.0,0.0))
    def check_mathieu_modcem2(self):
        cephes.mathieu_modcem2(1,1,1)
    def check_mathieu_sem(self):
        assert_equal(cephes.mathieu_sem(1,0,0),(0.0,1.0))
    def check_mathieu_modsem1(self):
        assert_equal(cephes.mathieu_modsem1(1,0,0),(0.0,0.0))
    def check_mathieu_modsem2(self):
        cephes.mathieu_modsem2(1,1,1)

    def check_modfresnelm(self):
        cephes.modfresnelm(0)
    def check_modfresnelp(self):
        cephes.modfresnelp(0)
    def _check_modstruve(self):
        assert_equal(cephes.modstruve(1,0),0.0)

    def check_nbdtr(self):
        assert_equal(cephes.nbdtr(1,1,1),1.0)
    def check_nbdtrc(self):
        assert_equal(cephes.nbdtrc(1,1,1),0.0)
    def check_nbdtri(self):
        assert_equal(cephes.nbdtri(1,1,1),1.0)
    def __check_nbdtrik(self):
        cephes.nbdtrik(1,.4,.5)
    def check_nbdtrin(self):
        assert_equal(cephes.nbdtrin(1,0,0),5.0)

    def check_ncfdtr(self):
        assert_equal(cephes.ncfdtr(1,1,1,0),0.0)
    def check_ncfdtri(self):
        assert_equal(cephes.ncfdtri(1,1,1,0),0.0)
    def check_ncfdtridfd(self):
        cephes.ncfdtridfd(1,0.5,0,1)
    def __check_ncfdtridfn(self):
        cephes.ncfdtridfn(1,0.5,0,1)
    def __check_ncfdtrinc(self):
        cephes.ncfdtrinc(1,0.5,0,1)

    def check_nctdtr(self):
        assert_equal(cephes.nctdtr(1,0,0),0.5)
    def __check_nctdtridf(self):
        cephes.nctdtridf(1,0.5,0)
    def check_nctdtrinc(self):
        cephes.nctdtrinc(1,0,0)
    def check_nctdtrit(self):
        cephes.nctdtrit(.1,0.2,.5)

    def check_ndtr(self):
        assert_equal(cephes.ndtr(0),0.5)
    def check_ndtri(self):
        assert_equal(cephes.ndtri(0.5),0.0)
    def check_nrdtrimn(self):
        assert_approx_equal(cephes.nrdtrimn(0.5,1,1),1.0)
    def check_nrdtrisd(self):
        assert_equal(cephes.nrdtrisd(0.5,0.5,0.5),0.0)

    def check_obl_ang1(self):
        cephes.obl_ang1(1,1,1,0)
    def check_obl_ang1_cv(self):
        result = cephes.obl_ang1_cv(1,1,1,1,0)
        assert_almost_equal(result[0],1.0)
        assert_almost_equal(result[1],0.0)

    def _check_obl_cv(self):
        assert_equal(cephes.obl_cv(1,1,0),2.0)
    def check_obl_rad1(self):
        cephes.obl_rad1(1,1,1,0)
    def check_obl_rad1_cv(self):
        cephes.obl_rad1_cv(1,1,1,1,0)
    def check_obl_rad2(self):
        cephes.obl_rad2(1,1,1,0)
    def check_obl_rad2_cv(self):
        cephes.obl_rad2_cv(1,1,1,1,0)

    def check_pbdv(self):
        assert_equal(cephes.pbdv(1,0),(0.0,0.0))
    def check_pbvv(self):
        cephes.pbvv(1,0)
    def check_pbwa(self):
        cephes.pbwa(1,0)
    def check_pdtr(self):
        cephes.pdtr(0,1)
    def check_pdtrc(self):
        cephes.pdtrc(0,1)
    def check_pdtri(self):
        cephes.pdtri(0.5,0.5)
    def check_pdtrik(self):
        cephes.pdtrik(0.5,1)

    def check_pro_ang1(self):
        cephes.pro_ang1(1,1,1,0)
    def check_pro_ang1_cv(self):
        assert_array_almost_equal(cephes.pro_ang1_cv(1,1,1,1,0),
                                  array((1.0,0.0)))
    def _check_pro_cv(self):
        assert_equal(cephes.pro_cv(1,1,0),2.0)
    def check_pro_rad1(self):
        cephes.pro_rad1(1,1,1,0.1) 
    def check_pro_rad1_cv(self):
        cephes.pro_rad1_cv(1,1,1,1,0)
    def check_pro_rad2(self):
        cephes.pro_rad2(1,1,1,0)
    def check_pro_rad2_cv(self):
        cephes.pro_rad2_cv(1,1,1,1,0)

    def check_psi(self):
        cephes.psi(1)

    def check_radian(self):
        assert_equal(cephes.radian(0,0,0),0)
    def check_rgamma(self):
        assert_equal(cephes.rgamma(1),1.0)
    def check_round(self):
        assert_equal(cephes.round(3.4),3.0)
        assert_equal(cephes.round(-3.4),-3.0)
        assert_equal(cephes.round(3.6),4.0)
        assert_equal(cephes.round(-3.6),-4.0)
        assert_equal(cephes.round(3.5),4.0)
        assert_equal(cephes.round(-3.5),-4.0)

    def check_shichi(self):
        cephes.shichi(1)
    def check_sici(self):
        cephes.sici(1)
    def check_sindg(self):
        assert_equal(cephes.sindg(90),1.0)
    def check_smirnov(self):
        assert_equal(cephes.smirnov(1,.1),0.9)
    def check_smirnovi(self):
        assert_almost_equal(cephes.smirnov(1,cephes.smirnovi(1,0.4)),0.4)
        assert_almost_equal(cephes.smirnov(1,cephes.smirnovi(1,0.6)),0.6)

    def check_spence(self):
        assert_equal(cephes.spence(1),0.0)
    def check_stdtr(self):
        assert_equal(cephes.stdtr(1,0),0.5)
    def check_stdtridf(self):
        cephes.stdtridf(0.7,1)
    def check_stdtrit(self):
        cephes.stdtrit(1,0.7)
    def check_struve(self):
        assert_equal(cephes.struve(0,0),0.0)

    def check_tandg(self):
        assert_equal(cephes.tandg(45),1.0)
    def check_tklmbda(self):
        assert_almost_equal(cephes.tklmbda(1,1),1.0)

    def check_y0(self):
        cephes.y0(1)
    def check_y1(self):
        cephes.y1(1)
    def check_yn(self):
        cephes.yn(1,1)
    def check_yv(self):
        cephes.yv(1,1)
    def _check_yve(self):
        cephes.yve(1,1)

    def check_zeta(self):
        cephes.zeta(2,2)
    def check_zetac(self):
        assert_equal(cephes.zetac(0),-1.5)
    def check_wofz(self):
        cephes.wofz(0)

class test_airy(unittest.TestCase):

    def check_airy(self):
	#This tests the airy function to ensure 8 place accuracy in computation

	x = airy(.99)
	assert_array_almost_equal(x,array([0.13689066,
					   -0.16050153,
					   1.19815925,
					   0.92046818]),8)
	x = airy(.41)
	assert_array_almost_equal(x,array([0.25238916,
					  -.23480512,
					  0.80686202,
					  0.51053919]),8)
	x = airy(-.36)
	assert_array_almost_equal(x,array([0.44508477,
					   -0.23186773,
					   0.44939534,
					   0.48105354]),8)

class test_airye(unittest.TestCase):

    def check_airye(self):
	a = airye(0.01)
	b = airy(0.01)
        b1 = [None]*4
        for n in range(2):
            b1[n] = b[n]*exp(2.0/3.0*0.01*sqrt(0.01))
        for n in range(2,4):
            b1[n] = b[n]*exp(-abs(real(2.0/3.0*0.01*sqrt(0.01))))
        assert_array_almost_equal(a,b1,6)

class test_arange(unittest.TestCase):

    def check_arange(self):
        numstring = arange(0,2.21,.1)
        assert_equal(numstring,array([0.,0.1,0.2,0.3,
                                      0.4,0.5,0.6,0.7,
                                      0.8,0.9,1.,1.1,
                                      1.2,1.3,1.4,1.5,
                                      1.6,1.7,1.8,1.9,
                                      2.,2.1,2.2]))
        numstringa = arange(3,4,.3)
        assert_array_equal(numstringa, array([3.,3.3,3.6,3.9]))
        numstringb = arange(3,27,3)
        assert_array_equal(numstringb,array([3,6,9,12,
                                             15,18,21,24]))
        numstringc = arange(3.3,27,4)
        assert_array_equal(numstringc,array([3.3,7.3,11.3,15.3,
                                             19.3,23.3]))

class test_ai_zeros(unittest.TestCase):

    def check_ai_zeros(self):
        ai = ai_zeros(1)
        assert_array_almost_equal(ai,(array([-2.33810741]),
                                     array([-1.01879297]),
                                     array([ 0.5357]),
                                     array([ 0.7012])),4)

class test_array(unittest.TestCase):

    def check_array(self):
        x = array([1,2,3,4])
        y = array([1,2,3,4])
        z = x*y
        assert_array_equal(z,array([1,4,9,16]))
        a = arange(1,5,1)
        assert_array_equal(a,x)

class test_assoc_laguerre(unittest.TestCase):

    def check_assoc_laguerre(self):
        a1 = genlaguerre(11,1)
        a2 = assoc_laguerre(.2,11,1)
        assert_array_almost_equal(a2,a1(.2),8)
        a2 = assoc_laguerre(1,11,1)
        assert_array_almost_equal(a2,a1(1),8)	

class test_besselpoly(unittest.TestCase):

    def check_besselpoly(self):
        pass

class test_bei(unittest.TestCase):

    def check_bei(self):
        mbei = bei(2)
        assert_almost_equal(mbei, 0.9722916273066613,5)#this may not be exact

class test_beip(unittest.TestCase):

    def check_beip(self):
        mbeip = beip(2)
        assert_almost_equal(mbeip,0.91701361338403631,5)#this may not be exact

class test_ber(unittest.TestCase):

    def check_ber(self):
        mber = ber(2)
        assert_almost_equal(mber,0.75173418271380821,5)#this may not be exact

class test_berp(unittest.TestCase):

    def check_berp(self):
        mberp = berp(2)
        assert_almost_equal(mberp,-0.49306712470943909,5)#this may not be exact

class test_bei_zeros(unittest.TestCase):

    def check_bei_zeros(self):
        bi = bi_zeros(5)
        assert_array_almost_equal(bi[0],array([-1.173713222709127,
                                               -3.271093302836352,
                                               -4.830737841662016,
                                               -6.169852128310251,
                                               -7.376762079367764]),11)

        assert_array_almost_equal(bi[1],array([-2.294439682614122,
                                               -4.073155089071828,
                                               -5.512395729663599,
                                               -6.781294445990305,
                                               -7.940178689168587]),10)
        
        assert_array_almost_equal(bi[2],array([-0.454944383639657,
                                               0.396522836094465,
                                               -0.367969161486959,
                                               0.349499116831805,
                                               -0.336026240133662]),11)
        
        assert_array_almost_equal(bi[3],array([0.601957887976239,
                                               -0.760310141492801,
                                               0.836991012619261,
                                               -0.88947990142654,
                                               0.929983638568022]),11)


class test_beip_zeros(unittest.TestCase):

    def check_beip_zeros(self):
        bip = beip_zeros(5)
        assert_array_almost_equal(bip,array([  3.772673304934953,
                                               8.280987849760042,
                                               12.742147523633703,
                                               17.193431752512542,
                                               21.641143941167325]),4)
class test_ber_zeros(unittest.TestCase):

    def check_ber_zeros(self):
        ber = ber_zeros(5)
        assert_array_almost_equal(ber,array([2.84892,
                                             7.23883,
                                             11.67396,
                                             16.11356,
                                             20.55463]),4)

class test_bernoulli(unittest.TestCase):

    def check_bernoulli(self):
        brn = bernoulli(5)
        assert_array_almost_equal(brn,array([1.0000,
                                             -0.5000,
                                             0.1667,
                                             0.0000,
                                             -0.0333,
                                             0.0000]),4)

class test_berp_zeros(unittest.TestCase):

    def check_berp_zeros(self):
        brp = berp_zeros(5)
        assert_array_almost_equal(brp,array([6.03871,
                                             10.51364,
                                             14.96844,
                                             19.41758,
                                             23.86430]),4)
class test_beta(unittest.TestCase):

    def check_beta(self):
        bet = beta(2,4)
        betg = (gamma(2)*gamma(4))/gamma(6)
        assert_almost_equal(bet,betg,8)

class test_betaln(unittest.TestCase):

    def check_betaln(self):
        betln = betaln(2,4)
        bet = log(abs(beta(2,4)))
        assert_almost_equal(betln,bet,8)

class test_betainc(unittest.TestCase):

    def check_betainc(self):
        btinc = betainc(1,1,.2)
        assert_almost_equal(btinc,0.2,8)

class test_betaincinv(unittest.TestCase):

    def check_betaincinv(self):
        y = betaincinv(2,4,.5)
        comp = betainc(2,4,y)
        assert_almost_equal(comp,.5,5)

class test_bi_zeros(unittest.TestCase):

    def check_bi_zeros(self):
        bi = bi_zeros(2)
        bia = (array([-1.17371322, -3.2710930]),
        array([-2.29443968, -4.07315509]),
        array([-0.45494438,  0.39652284]),
        array([ 0.60195789 , -0.76031014]))
        assert_array_almost_equal(bi,bia,4)

class test_chebyc(unittest.TestCase):

    def check_chebyc(self):
        C0 = chebyc(0)
        C1 = chebyc(1)
        C2 = chebyc(2)
        C3 = chebyc(3)
        C4 = chebyc(4)
        C5 = chebyc(5)

        assert_array_almost_equal(C0.c,[2],13)
        assert_array_almost_equal(C1.c,[1,0],13)
        assert_array_almost_equal(C2.c,[1,0,-2],13)
        assert_array_almost_equal(C3.c,[1,0,-3,0],13)
        assert_array_almost_equal(C4.c,[1,0,-4,0,2],13)
        assert_array_almost_equal(C5.c,[1,0,-5,0,5,0],13)

class test_chebys(unittest.TestCase):

    def check_chebys(self):
        S0 = chebys(0)
        S1 = chebys(1)        
        S2 = chebys(2)
        S3 = chebys(3)
        S4 = chebys(4)
        S5 = chebys(5)
        assert_array_almost_equal(S0.c,[1],13)
        assert_array_almost_equal(S1.c,[1,0],13)
        assert_array_almost_equal(S2.c,[1,0,-1],13)
        assert_array_almost_equal(S3.c,[1,0,-2,0],13)
        assert_array_almost_equal(S4.c,[1,0,-3,0,1],13)
        assert_array_almost_equal(S5.c,[1,0,-4,0,3,0],13)

class test_chebyt(unittest.TestCase):

    def check_chebyt(self):
        T0 = chebyt(0)
        T1 = chebyt(1)        
        T2 = chebyt(2)
        T3 = chebyt(3)
        T4 = chebyt(4)
        T5 = chebyt(5)
        assert_array_almost_equal(T0.c,[1],13)
        assert_array_almost_equal(T1.c,[1,0],13)
        assert_array_almost_equal(T2.c,[2,0,-1],13)
        assert_array_almost_equal(T3.c,[4,0,-3,0],13)
        assert_array_almost_equal(T4.c,[8,0,-8,0,1],13)
        assert_array_almost_equal(T5.c,[16,0,-20,0,5,0],13)

class test_chebyu(unittest.TestCase):

    def check_chebyu(self):
        U0 = chebyu(0)
        U1 = chebyu(1)        
        U2 = chebyu(2)
        U3 = chebyu(3)
        U4 = chebyu(4)
        U5 = chebyu(5)
        assert_array_almost_equal(U0.c,[1],13)
        assert_array_almost_equal(U1.c,[2,0],13)
        assert_array_almost_equal(U2.c,[4,0,-1],13)
        assert_array_almost_equal(U3.c,[8,0,-4,0],13)
        assert_array_almost_equal(U4.c,[16,0,-12,0,1],13)
        assert_array_almost_equal(U5.c,[32,0,-32,0,6,0],13)

class test_choose(unittest.TestCase):

    def check_choose(self):
        carray = [1,3,2,4,6,5]
        chose = choose([1,3,5],carray)
        assert_array_equal(chose,array([3,4,5]))

class test_cbrt(unittest.TestCase):

    def check_cbrt(self):
        cb = cbrt(27)
        cbrl = 27**(1.0/3.0)
        assert_approx_equal(cb,cbrl)

    def check_cbrtmore(self):
        cb1 = cbrt(27.9)
        cbrl1 = 27.9**(1.0/3.0)
        assert_almost_equal(cb1,cbrl1,8)

class test_cosdg(unittest.TestCase):

    def check_cosdg(self):
        cdg = cosdg(90)
        cdgrl = cos(pi/2.0)
        assert_almost_equal(cdg,cdgrl,8)

    def check_cosdgmore(self):
        cdgm = cosdg(30)
        cdgmrl = cos(pi/6.0)
        assert_almost_equal(cdgm,cdgmrl,8)

class test_cosm1(unittest.TestCase):

    def check_cosm1(self):
        cs = (cosm1(0),cosm1(.3),cosm1(pi/10))
        csrl = (cos(0)-1,cos(.3)-1,cos(pi/10)-1)
        assert_array_almost_equal(cs,csrl,8)

class test_cotdg(unittest.TestCase):

    def check_cotdg(self):
        ct = cotdg(30)
        ctrl = tan(pi/6.0)**(-1)
        assert_almost_equal(ct,ctrl,8)

    def check_cotdgmore(self):
        ct1 = cotdg(45)
        ctrl1 = tan(pi/4.0)**(-1)
        assert_almost_equal(ct1,ctrl1,8)

class test_ellipj(unittest.TestCase):

    def check_ellipj(self):
        el = ellipj(0.2,0)
        rel = [sin(0.2),cos(0.2),1.0,0.20]
        assert_array_almost_equal(el,rel,13)

class test_ellipk(unittest.TestCase):

    def check_ellipk(self):
        elk = ellipk(.2)
        assert_almost_equal(elk,1.659623598610528,11)

class test_ellipkinc(unittest.TestCase):

    def check_ellipkinc(self):
        elkinc = ellipkinc(pi/2,.2)
        elk = ellipk(0.2)
        assert_almost_equal(elkinc,elk,15)
        alpha = 20*pi/180
        phi = 45*pi/180
        m = sin(alpha)**2
        elkinc = ellipkinc(phi,m)
        assert_almost_equal(elkinc,0.79398143,8)
        # From pg. 614 of A & S


class test_ellipe(unittest.TestCase):

    def check_ellipe(self):
        ele = ellipe(.2)
        assert_almost_equal(ele,1.4890350580958529,8)

class test_ellipeinc(unittest.TestCase):

    def check_ellipeinc(self):
        eleinc = ellipeinc(pi/2,.2)
        ele = ellipe(0.2)
        assert_almost_equal(eleinc,ele,14)
        # pg 617 of A & S
        alpha, phi = 52*pi/180,35*pi/180
        m = sin(alpha)**2
        eleinc = ellipeinc(phi,m)
        assert_almost_equal(eleinc, 0.58823065, 8)
        

class test_erf(unittest.TestCase):

    def check_erf(self):
        er = erf(.25)
        assert_almost_equal(er,0.2763263902,8)

class test_erf_zeros(unittest.TestCase):

    def check_erf_zeros(self):
        erz = erf_zeros(5)
        erzr= array([1.45061616+1.88094300j,
                     2.24465928+2.61657514j,
                     2.83974105+3.17562810j,
                     3.33546074+3.64617438j,
                     3.76900557+4.06069723j])
        assert_array_almost_equal(erz,erzr,4)

class test_erfcinv(unittest.TestCase):

    def check_erfcinv(self):
        i = erfcinv(1)
        assert_equal(i,0)

class test_erfinv(unittest.TestCase):

    def check_erfinv(self):
        i = erfinv(0)
        assert_equal(i,0)

class test_errprint(unittest.TestCase):

    def check_errprint(self):
        a = errprint()
        b = 1-a #a is the state 1-a inverts state
        c = errprint(b) #returns last state 'a' 
        assert_equal(a,c)
        d = errprint(a) #returns to original state
        assert_equal(d,b) #makes sure state was returned
        #assert_equal(d,1-a)

class test_euler(unittest.TestCase):

    def check_euler(self):
        eu0 = euler(0)
        eu1 = euler(1)  
        eu2 = euler(2)   # just checking segfaults
        assert_almost_equal(eu0[0],1,8)
        assert_almost_equal(eu2[2],-1,8)
        eu24 = euler(24)
        mathworld = [1,1,5,61,1385,50521,2702765,199360981,
                     19391512145l,2404879675441l,
                     370371188237525l,69348874393137901l,
                     15514534163557086905l]
        correct = zeros((25,),'d')
        for k in range(0,13):
            if (k % 2):
                correct[2*k] = -float(mathworld[k])
            else:
                correct[2*k] = float(mathworld[k])
        err = nan_to_num((eu24-correct)/correct)
        errmax = max(err)
        assert_almost_equal(errmax, 0.0, 14)

class test_exp2(unittest.TestCase):

    def check_exp2(self):
        ex = exp2(2)
        exrl = 2**2
        assert_equal(ex,exrl)

    def check_exp2more(self):
        exm = exp2(2.5)
        exmrl = 2**(2.5)
        assert_almost_equal(exm,exmrl,8)

class test_exp10(unittest.TestCase):

    def check_exp10(self):
        ex = exp10(2)
        exrl = 10**2
        assert_approx_equal(ex,exrl)

    def check_exp10more(self):
        exm = exp10(2.5)
        exmrl = 10**(2.5)
        assert_almost_equal(exm,exmrl,8)

class test_expm1(unittest.TestCase):

    def check_expm1(self):
        ex = (expm1(2),expm1(3),expm1(4))
        exrl = (exp(2)-1,exp(3)-1,exp(4)-1)
        assert_array_almost_equal(ex,exrl,8)

    def check_expm1more(self):
        ex1 = (expm1(2),expm1(2.1),expm1(2.2))
        exrl1 = (exp(2)-1,exp(2.1)-1,exp(2.2)-1)
        assert_array_almost_equal(ex1,exrl1,8)

class test_fresnel(unittest.TestCase):

    def check_fresnel(self):
        frs = array(fresnel(.5))
        assert_array_almost_equal(frs,array([0.064732432859999287, 0.49234422587144644]),8)

class test_fresnel_zeros(unittest.TestCase):

    # values from pg 329  Table 7.11 of A & S
    #  slightly corrected in 4th decimal place 
    def check_fresnel_zeros(self):
        szo, czo = fresnel_zeros(5)
        assert_array_almost_equal(szo,
                                  array([ 2.0093+0.2885j,
                                          2.8335+0.2443j,
                                          3.4675+0.2185j,
                                          4.0026+0.2009j,
                                          4.4742+0.1877j]),3)
        assert_array_almost_equal(czo,
                                  array([ 1.7437+0.3057j,
                                          2.6515+0.2529j,
                                          3.3204+0.2240j,
                                          3.8757+0.2047j,
                                          4.3611+0.1907j]),3)
        vals1 = fresnel(szo)[0]
        vals2 = fresnel(czo)[1]
        assert_array_almost_equal(vals1,0,14)
        assert_array_almost_equal(vals2,0,14)
        
    def check_fresnelc_zeros(self):
        szo, czo = fresnel_zeros(6)
        frc = fresnelc_zeros(6)
        assert_array_almost_equal(frc,czo,12)

    def check_fresnels_zeros(self):
        szo, czo = fresnel_zeros(5)
        frs = fresnels_zeros(5)
        assert_array_almost_equal(frs,szo,12)


class test_gamma(unittest.TestCase):

    def check_gamma(self):

        gam = gamma(5)
        assert_equal(gam,24.0)

class test_gammaln(unittest.TestCase):

    def check_gammaln(self):
        gamln = gammaln(3)
        lngam = log(gamma(3))
        assert_almost_equal(gamln,lngam,8)

class test_gammainc(unittest.TestCase):

    def check_gammainc(self):
        gama = gammainc(.5,.5)
        assert_almost_equal(gama,.7,1)

class test_gammaincc(unittest.TestCase):

    def check_gammaincc(self):
        gicc = gammaincc(.5,.5)
        greal = 1 - gammainc(.5,.5)
        assert_almost_equal(gicc,greal,8)

class test_gammainccinv(unittest.TestCase):

    def check_gammainccinv(self):
        gccinv = gammainccinv(.5,.5)
        gcinv = gammaincinv(.5,.5)
        assert_almost_equal(gccinv,gcinv,8)

class test_gammaincinv(unittest.TestCase):

    def check_gammaincinv(self):
        y = gammaincinv(.4,.4)
        x = gammainc(.4,y)
        assert_almost_equal(x,0.4,1)

class test_hankel1(unittest.TestCase):

    def check_hankel1(self):
        hank1 = hankel1(1,.1)
        hankrl = (jv(1,.1)+yv(1,.1)*1j)
        assert_almost_equal(hank1,hankrl,8)

class test_hankel1e(unittest.TestCase):

    def check_hankel1e(self):
        hank1e = hankel1e(1,.1)
        hankrle = hankel1(1,.1)*exp(-.1j)
        assert_almost_equal(hank1e,hankrle,8)

class test_hankel2(unittest.TestCase):

    def check_hankel2(self):
        hank2 = hankel2(1,.1)
        hankrl2 = (jv(1,.1)-yv(1,.1)*1j)
        assert_almost_equal(hank2,hankrl2,8)

class test_hankel2e(unittest.TestCase):

    def check_hankl2e(self):
        hank2e = hankel2e(1,.1)
        hankrl2e = hankel2e(1,.1)
        assert_almost_equal(hank2e,hankrl2e,8)

class test_hermite(unittest.TestCase):

    def check_hermite(self):
        H0 = hermite(0)
        H1 = hermite(1)        
        H2 = hermite(2)
        H3 = hermite(3)
        H4 = hermite(4)
        H5 = hermite(5)
        assert_array_almost_equal(H0.c,[1],13)
        assert_array_almost_equal(H1.c,[2,0],13)
        assert_array_almost_equal(H2.c,[4,0,-2],13)
        assert_array_almost_equal(H3.c,[8,0,-12,0],13)
        assert_array_almost_equal(H4.c,[16,0,-48,0,12],12)
        assert_array_almost_equal(H5.c,[32,0,-160,0,120,0],12)

    def check_hermitenorm(self):
        # He_n(x) = 2**(-n/2) H_n(x/sqrt(2))
        psub = poly1d([1.0/sqrt(2),0])
        H0 = hermitenorm(0)
        H1 = hermitenorm(1)        
        H2 = hermitenorm(2)
        H3 = hermitenorm(3)
        H4 = hermitenorm(4)
        H5 = hermitenorm(5)
        he0 = hermite(0)(psub)
        he1 = hermite(1)(psub) / sqrt(2)
        he2 = hermite(2)(psub) / 2.0
        he3 = hermite(3)(psub) / (2*sqrt(2))
        he4 = hermite(4)(psub) / 4.0
        he5 = hermite(5)(psub) / (4.0*sqrt(2))

        assert_array_almost_equal(H0.c,he0.c,13)
        assert_array_almost_equal(H1.c,he1.c,13)
        assert_array_almost_equal(H2.c,he2.c,13)
        assert_array_almost_equal(H3.c,he3.c,13)
        assert_array_almost_equal(H4.c,he4.c,13)
        assert_array_almost_equal(H5.c,he5.c,13)

_gam = cephes.gamma
    
class test_gegenbauer(unittest.TestCase):

    def check_gegenbauer(self):
        a = 5*rand()-0.5
        if (a==0): a = -0.2
        print "Gegenbauer, a = ", a
        Ca0 = gegenbauer(0,a)
        Ca1 = gegenbauer(1,a)
        Ca2 = gegenbauer(2,a)
        Ca3 = gegenbauer(3,a)
        Ca4 = gegenbauer(4,a)
        Ca5 = gegenbauer(5,a)

        assert_array_almost_equal(Ca0.c,array([1]),13)
        assert_array_almost_equal(Ca1.c,array([2*a,0]),13)
        assert_array_almost_equal(Ca2.c,array([2*a*(a+1),0,-a]),13)
        assert_array_almost_equal(Ca3.c,array([4*poch(a,3),0,-6*a*(a+1),
                                               0])/3.0,11)
        assert_array_almost_equal(Ca4.c,array([4*poch(a,4),0,-12*poch(a,3),
                                               0,3*a*(a+1)])/6.0,11)
        assert_array_almost_equal(Ca5.c,array([4*poch(a,5),0,-20*poch(a,4),
                                               0,15*poch(a,3),0])/15.0,11)
        

class test_h1vp(unittest.TestCase):

    def check_h1vp(self):
        
        h1 = h1vp(1,.1)
        h1real = (jvp(1,.1)+yvp(1,.1)*1j)
        assert_almost_equal(h1,h1real,8)

class test_h2vp(unittest.TestCase):

    def check_h2vp(self):
        h2 = h2vp(1,.1)
        h2real = (jvp(1,.1)-yvp(1,.1)*1j)
        assert_almost_equal(h2,h2real,8)

class test_hyp0f1(unittest.TestCase):

    def check_hyp0f1(self):
        pass

class test_hyp1f1(unittest.TestCase):

    def check_hyp1f1(self):

        hyp1 = hyp1f1(.1,.1,.3)
        assert_almost_equal(hyp1, 1.3498588075760032,7)

class test_hyp1f2(unittest.TestCase):

    def check_hyp1f2(self):
        pass

class test_hyp2f0(unittest.TestCase):

    def check_hyp2f0(self):
        pass

class test_hyp2f1(unittest.TestCase):

    def check_hyp2f1(self):
        hyp = hyp2f1(1,1,2,.5)
        hrl = -(1/.5)*log(1-.5)
        assert_almost_equal(hyp,hrl,8)

class test_hyp3f0(unittest.TestCase):

    def check_hyp3f0(self):
        pass

class test_hyperu(unittest.TestCase):

    def check_hyperu(self):
        val1 = hyperu(1,0.1,100)
        assert_almost_equal(val1,0.0098153,7)
        a,b = [0.3,0.6,1.2,-2.7],[1.5,3.2,-0.4,-3.2]
        a,b = asarray(a), asarray(b)
        z = 0.5
        hypu = hyperu(a,b,z)
        hprl = (pi/sin(pi*b))*(hyp1f1(a,b,z)/ \
                               (gamma(1+a-b)*gamma(b))- \
                               z**(1-b)*hyp1f1(1+a-b,2-b,z) \
                               /(gamma(a)*gamma(2-b)))
        assert_array_almost_equal(hypu,hprl,12)

class test_i0(unittest.TestCase):

    def check_i0(self):
        oiz = i0(.1)
        oizr = iv(0,.1)
        assert_almost_equal(oiz,oizr,8)

class test_i0e(unittest.TestCase):

    def check_i0e(self):
        oize = i0e(.1)
        oizer = ive(0,.1)
        assert_almost_equal(oize,oizer,8)

class test_i1(unittest.TestCase):

    def check_i1(self):

        oi1 = i1(.1)
        oi1r = iv(1,.1)
        assert_almost_equal(oi1,oi1r,8)

class test_i1e(unittest.TestCase):

    def check_i1e(self):
        oi1e = i1e(.1)
        oi1er = ive(1,.1)
        assert_almost_equal(oi1e,oi1er,8)

class test_iti0k0(unittest.TestCase):

    def check_iti0k0(self):
        iti0 = array(iti0k0(5))
        assert_array_almost_equal(iti0,array([31.848667776169801, 1.5673873907283657]),5)

class test_it2i0k0(unittest.TestCase):

    def check_it2i0k0(self):
        it2k = it2i0k0(.1)
        assert_array_almost_equal(it2k,array([0.0012503906973464409, 3.3309450354686687]),6)

class test_itj0y0(unittest.TestCase):

    def check_itj0y0(self):
        it0 = array(itj0y0(.2))
        assert_array_almost_equal(it0,array([0.19933433254006822, -0.34570883800412566]),8)

class test_it2j0y0(unittest.TestCase):

    def check_it2j0y0(self):
        it2 = array(it2j0y0(.2))
        assert_array_almost_equal(it2,array([0.0049937546274601858, -0.43423067011231614]),8)

class test_iv(unittest.TestCase):

    def check_iv(self):
        iv1 = iv(0,.1)*exp(-.1)
        assert_almost_equal(iv1,0.90710092578230106,10)

class test_ive(unittest.TestCase):

    def check_ive(self):
        ive1 = ive(0,.1)
        iv1 = iv(0,.1)*exp(-.1)
        assert_almost_equal(ive1,iv1,10)

class test_ivp(unittest.TestCase):

    def check_ivp(self):
        y=(iv(0,2)-iv(2,2))/2
        x = ivp(1,2)
        assert_almost_equal(x,y,10)

class test_j0(unittest.TestCase):

    def check_j0(self):
        oz = j0(.1)
        ozr = jn(0,.1)
        assert_almost_equal(oz,ozr,8)

class test_j1(unittest.TestCase):

    def check_j1(self):
        o1 = j1(.1)
        o1r = jn(1,.1)
        assert_almost_equal(o1,o1r,8)

class test_jacobi(unittest.TestCase):

    def check_jacobi(self):
        a = 5*rand() - 1
        b = 5*rand() - 1
        P0 = jacobi(0,a,b)
        P1 = jacobi(1,a,b)
        P2 = jacobi(2,a,b)
        P3 = jacobi(3,a,b)
        
        assert_array_almost_equal(P0.c,[1],13)
        assert_array_almost_equal(P1.c,array([a+b+2,a-b])/2.0,13)
        cp = [(a+b+3)*(a+b+4), 4*(a+b+3)*(a+2), 4*(a+1)*(a+2)]
        p2c = [cp[0],cp[1]-2*cp[0],cp[2]-cp[1]+cp[0]]
        assert_array_almost_equal(P2.c,array(p2c)/8.0,13)
        cp = [(a+b+4)*(a+b+5)*(a+b+6),6*(a+b+4)*(a+b+5)*(a+3),
              12*(a+b+4)*(a+2)*(a+3),8*(a+1)*(a+2)*(a+3)]
        p3c = [cp[0],cp[1]-3*cp[0],cp[2]-2*cp[1]+3*cp[0],cp[3]-cp[2]+cp[1]-cp[0]]
        assert_array_almost_equal(P3.c,array(p3c)/48.0,13)        
        

class test_jn(unittest.TestCase):

    def check_jn(self):
        jnnr = jn(1,.2)
        assert_almost_equal(jnnr,0.099500832639235995,8)

class test_jv(unittest.TestCase):

    def check_jv(self):
        jc = jv(0,.1)
        assert_almost_equal(jc,0.99750156206604002,8)

class test_jve(unittest.TestCase):

    def check_jve(self):
        jvexp = jve(1,.2)
        assert_almost_equal(jvexp,0.099500832639235995,8)
        jvexp1 = jve(1,.2+1j)
        z = .2+1j
        jvexpr = jv(1,z)*exp(-abs(z.imag))
        assert_almost_equal(jvexp1,jvexpr,8)

class test_jn_zeros(unittest.TestCase):

    def check_jn_zeros(self):
        jn0 = jn_zeros(0,5)
        jn1 = jn_zeros(1,5)
        assert_array_almost_equal(jn0,array([ 2.4048255577,
                                              5.5200781103,
                                              8.6537279129,
                                              11.7915344391,
                                              14.9309177086]),4)
        assert_array_almost_equal(jn1,array([ 3.83171,
                                              7.01559,
                                              10.17347,
                                              13.32369,
                                              16.47063]),4)

class test_jnjnp_zeros(unittest.TestCase):

    def check_jnjnp_zeros(self):
        pass
        #jnjp = jnjnp(3)
        #assert_array_almost_equal(jnjp,(array([
        #I don't think specfun jdzo is working properly the outputs do not seem to correlate
        #to the inputs

class test_jnp_zeros(unittest.TestCase):

    def check_jnp_zeros(self):
        jnp = jnp_zeros(1,5)
        assert_array_almost_equal(jnp, array([	1.84118,
                                                5.33144,
                                                8.53632,
                                                11.70600,
                                                14.86359]),4)

class test_jnyn_zeros(unittest.TestCase):

    def check_jnyn_zeros(self):
        jnz = jnyn_zeros(1,5)
        assert_array_almost_equal(jnz,(array([	3.83171,
                                                7.01559,
                                                10.17347,
                                                13.32369,
                                                16.47063]),
                                       array([	1.84118,
                                                5.33144,
                                                8.53632,
                                                11.70600,
                                                14.86359]),
                                       array([	2.19714,
                                                5.42968,
                                                8.59601,
                                                11.74915,
                                                14.89744]),
                                       array([	3.68302,
                                                6.94150,
                                                10.12340,
                                                13.28576,
                                                16.44006])),4)

class test_jvp(unittest.TestCase):

    def check_jvp(self):
        jvprim = jvp(2,2)
        jv0 = (jv(1,2)-jv(3,2))/2
        assert_almost_equal(jvprim,jv0,4)

class test_k0(unittest.TestCase):

    def check_k0(self):
        ozk = k0(.1)
        ozkr = kv(0,.1)
        assert_almost_equal(ozk,ozkr,8)

class test_k0e(unittest.TestCase):

    def check_k0e(self):
        ozke = k0e(.1)
        ozker = kve(0,.1)
        assert_almost_equal(ozke,ozker,8)

class test_k1(unittest.TestCase):

    def check_k1(self):
        o1k = k1(.1)
        o1kr = kv(1,.1)
        assert_almost_equal(o1k,o1kr,8)

class test_k1e(unittest.TestCase):

    def check_k1e(self):
        o1ke = k1e(.1)
        o1ker = kve(1,.1)
        assert_almost_equal(o1ke,o1ker,8)

class test_kei(unittest.TestCase):

    def check_kei(self):
        mkei = kei(2)
        assert_almost_equal(mkei,-0.20240006776470432,5)

class test_kelvin(unittest.TestCase):

    def check_kelvin(self):
        mkelv = kelvin(2)
        assert_array_almost_equal(mkelv,(ber(2)+bei(2)*1j,
                                         ker(2)+kei(2)*1j,
                                         berp(2)+beip(2)*1j,
                                         kerp(2)+keip(2)*1j),8)

class test_keip(unittest.TestCase):

    def  check_keip(self):
        mkeip = keip(2)
        assert_almost_equal(mkeip,0.21980790991960536,5)

class test_ker(unittest.TestCase):

    def check_ker(self):
        mker = ker(2)
        assert_almost_equal(mker,-0.041664513991509472,5)

class test_kerp(unittest.TestCase):

    def check_kerp(self):
        mkerp = kerp(2)
        assert_almost_equal(mkerp,-0.10660096588105264,5)

class test_kei_zeros(unittest.TestCase):

    def check_kei_zeros(self):
        kei = kei_zeros(5)
        assert_array_almost_equal(kei,array([  3.91467,
                                              8.34422,
                                              12.78256,
                                              17.22314,
                                              21.66464]),4)

class test_keip_zeros(unittest.TestCase):

    def check_keip_zeros(self):
        keip = keip_zeros(5)
        assert_array_almost_equal(keip,array([	4.93181,
                                                9.40405,
                                                13.85827,
                                                18.30717,
                                                22.75379]),4)



class test_kelvin_zeros(unittest.TestCase):

    # numbers come from 9.9 of A&S pg. 381 
    def check_kelvin_zeros(self):
        tmp = kelvin_zeros(5)
        berz,beiz,kerz,keiz,berpz,beipz,kerpz,keipz = tmp
        assert_array_almost_equal(berz,array([ 2.84892,
                                               7.23883,
                                               11.67396,
                                               16.11356,
                                               20.55463]),4)
        assert_array_almost_equal(beiz,array([ 5.02622,
                                               9.45541,
                                               13.89349,
                                               18.33398,
                                               22.77544]),4)
        assert_array_almost_equal(kerz,array([ 1.71854,
                                               6.12728,
                                               10.56294,
                                               15.00269,
                                               19.44382]),4)
        assert_array_almost_equal(keiz,array([ 3.91467,
                                               8.34422,
                                               12.78256,
                                               17.22314,
                                               21.66464]),4)
        assert_array_almost_equal(berpz,array([ 6.03871,
                                                10.51364,
                                                14.96844,
                                                19.41758,
                                                23.86430]),4)
        assert_array_almost_equal(beipz,array([ 3.77267,
                 # table from 1927 had 3.77320
                 #  but this is more accurate
                                                8.28099,
                                                12.74215,
                                                17.19343,
                                                21.64114]),4)
        assert_array_almost_equal(kerpz,array([	2.66584,
                                                7.17212,
                                                11.63218,
                                                16.08312,
                                                20.53068]),4)
        assert_array_almost_equal(keipz,array([	4.93181,
                                                9.40405,
                                                13.85827,
                                                18.30717,
                                                22.75379]),4)
        
class test_ker_zeros(unittest.TestCase):

    def check_ker_zeros(self):
        ker = ker_zeros(5)
        assert_array_almost_equal(ker,array([  1.71854,
                                               6.12728,
                                               10.56294,
                                               15.00269,
                                               19.44381]),4)

class test_kerp_zeros(unittest.TestCase):

    def check_kerp_zeros(self):
        kerp = kerp_zeros(5)
        assert_array_almost_equal(kerp,array([  2.66584,
                                                7.17212,
                                                11.63218,
                                                16.08312,
                                                20.53068]),4)

class test_kn(unittest.TestCase):

    def check_kn(self):
        kn1 = kn(0,.2)
        assert_almost_equal(kn1,1.7527038555281462,8)

class test_kv(unittest.TestCase):

    def check_kv(self):
        kv1 = kv(0,.2)
        assert_almost_equal(kv1,1.7527038555281462,8)

class test_kve(unittest.TestCase):

    def check_kve(self):
        kve1 = kve(0,.2)
        kv1 = kv(0,.2)*exp(.2)
        assert_almost_equal(kve1,kv1,8)
        z = .2+1j
        kve2 = kve(0,z)
        kv2 = kv(0,z)*exp(z)
        assert_almost_equal(kve2,kv2,8)

class test_kvp(unittest.TestCase):

    def check_kvp(self):
        kvprim = kvp(1,2)
        kvprimrl = (kv(0,2) - kv(2,2))/2
        assert_almost_equal(kvprim,kvprimrl,4)	 #this function (kvp) is broken

class test_laguerre(unittest.TestCase):

    def check_laguerre(self):
        lag0 = laguerre(0)
        lag1 = laguerre(1)
        lag2 = laguerre(2)
        lag3 = laguerre(3)
        lag4 = laguerre(4)
        lag5 = laguerre(5)
        assert_array_almost_equal(lag0.c,[1],13)
        assert_array_almost_equal(lag1.c,[-1,1],13)
        assert_array_almost_equal(lag2.c,array([1,-4,2])/2.0,13)
        assert_array_almost_equal(lag3.c,array([-1,9,-18,6])/6.0,13)
        assert_array_almost_equal(lag4.c,array([1,-16,72,-96,24])/24.0,13)
        assert_array_almost_equal(lag5.c,array([-1,25,-200,600,-600,120])/120.0,13)

    def check_genlaguerre(self):
        k = 5*rand()-0.9
        lag0 = genlaguerre(0,k)
        lag1 = genlaguerre(1,k)
        lag2 = genlaguerre(2,k)
        lag3 = genlaguerre(3,k)
        assert_equal(lag0.c,[1])
        assert_equal(lag1.c,[-1,k+1])
        assert_equal(lag2.c,array([1,-2*(k+2),(k+1.)*(k+2.)])/2.0)
        assert_equal(lag3.c,array([-1,3*(k+3),-3*(k+2)*(k+3),(k+1)*(k+2)*(k+3)])/6.0)
        

# Base polynomials come from Abrahmowitz and Stegan
class test_legendre(unittest.TestCase):

    def check_legendre(self):
        leg0 = legendre(0)
        leg1 = legendre(1)
        leg2 = legendre(2)
        leg3 = legendre(3)
        leg4 = legendre(4)
        leg5 = legendre(5)
        assert_equal(leg0.c,[1])
        assert_equal(leg1.c,[1,0])
        assert_equal(leg2.c,array([3,0,-1])/2.0)
        assert_equal(leg3.c,array([5,0,-3,0])/2.0)
        assert_equal(leg4.c,array([35,0,-30,0,3])/8.0)
        assert_equal(leg5.c,array([63,0,-70,0,15,0])/8.0)


class test_lmbda(unittest.TestCase):

    def check_lmbda(self):
        lam = lmbda(1,.1)
        lamr = (array([jn(0,.1), 2*jn(1,.1)/.1]),
                array([jvp(0,.1), -2*jv(1,.1)/.01 + 2*jvp(1,.1)/.1]))
        assert_array_almost_equal(lam,lamr,8)

class test_log1p(unittest.TestCase):

    def check_log1p(self):
        l1p = (log1p(10),log1p(11),log1p(12))
        l1prl = (log(11),log(12),log(13))
        assert_array_almost_equal(l1p,l1prl,8)

    def check_log1pmore(self):
        l1pm = (log1p(1),log1p(1.1),log1p(1.2))
        l1pmrl = (log(2),log(2.1),log(2.2))
        assert_array_almost_equal(l1pm,l1pmrl,8)

class test_lpmn(unittest.TestCase):

    def check_lpmn(self):

        lp = lpmn(0,2,.5)
        assert_array_almost_equal(lp,(array([	    [ 1.00000 ,
                                                      0.50000,
                                                      -0.12500]]),
                                      array([	    [ 0.00000 ,
                                                      1.00000 ,
                                                      1.50000]])),4)

class test_lpn(unittest.TestCase):

    def check_lpn(self):
        lpnf = lpn(2,.5)
        assert_array_almost_equal(lpnf,(array(      [ 1.00000 ,
                                                        0.50000,
                                                        -0.12500]),
                                      array(	    [ 0.00000 ,
                                                      1.00000 ,
                                                      1.50000])),4)

class test_lpmv(unittest.TestCase):

    def check_lpmv(self):
        lp = lpmv(0,2,.5)
        assert_almost_equal(lp,-0.125,3)

class test_lqmn(unittest.TestCase):

    def check_lqmn(self):
        lqmnf = lqmn(0,2,.5)
        lqmnf = lqmn(0,2,.5)
        lqf = lqn(2,.5)
        assert_array_almost_equal(lqmnf[0][0],lqf[0],4)
        assert_array_almost_equal(lqmnf[1][0],lqf[1],4)

                


class test_lqn(unittest.TestCase):

    def check_lqn(self):
        lqf = lqn(2,.5)
        assert_array_almost_equal(lqf,(array([ 0.5493, -0.7253, -0.8187]),
                                       array([ 1.3333,	1.216 , -0.8427])),4)

class test_mathieu_a(unittest.TestCase):

    def check_mathieu_a(self):
        pass

class test_mathieu_even_coef(unittest.TestCase):

    def check_mathieu_even_coef(self):
        mc =  mathieu_even_coef(2,5)
        #Q not defined broken and cannot figure out proper reporting order

class test_mathieu_odd_coef(unittest.TestCase):

    def check_mathieu_odd_coef(self):
        pass
            #same problem as above

class test_modfresnelp(unittest.TestCase):

    def check_modfresnelp(self):
        pass

class test_modfresnelm(unittest.TestCase):

    def check_modfresnelm(self):
        pass

class test_obl_cv_seq(unittest.TestCase):

    def check_obl_cv_seq(self):
        obl = obl_cv_seq(0,3,1)
        assert_array_almost_equal(obl,array([ -0.348602,
                                              1.393206,
                                              5.486800,
                                              11.492120]),5)

class test_pbdn_seq(unittest.TestCase):

    def check_pbdn_seq(self):
        pb = pbdn_seq(1,.1)
        assert_array_almost_equal(pb,(array([ 0.9975,
                                              0.0998]),
                                      array([-0.0499,
                                             0.9925])),4)

class test_pbdv(unittest.TestCase):

    def check_pbdv(self):
        pbv = pbdv(1,.2)
        derrl = 1/2*(.2)*pbdv(1,.2)[0] - pbdv(0,.2)[0]

class _test_pbdv_seq(unittest.TestCase):

    def check_pbdv_seq(self):
        pbn = pbdn_seq(1,.1)
        pbv = pbdv_seq(1,.1)
        assert_array_almost_equal(pbv,(real(pbn[0]),real(pbn[1])),4)

class test_pbvv_seq(unittest.TestCase):

    def check_pbvv_seq(self):
        pass

class test_polygamma(unittest.TestCase):

    # from Table 6.2 (pg. 271) of A&S
    def check_polygamma(self):
        poly2 = polygamma(2,1)
        poly3 = polygamma(3,1)
        assert_almost_equal(poly2,-2.4041138063,10)
        assert_almost_equal(poly3,6.4939394023,10)

class test_pro_cv_seq(unittest.TestCase):

    def check_pro_cv_seq(self):
        prol = pro_cv_seq(0,3,1)
        assert_array_almost_equal(prol,array([  0.319000,
                                               2.593084,
                                               6.533471,
                                               12.514462]),5)

class test_psi(unittest.TestCase):

    def check_psi(self):
        ps = psi(1)
        assert_almost_equal(ps,-0.57721566490153287,8)

class test_radian(unittest.TestCase):

    def check_radian(self):
        rad = radian(90,0,0)
        assert_almost_equal(rad,pi/2.0,5)

    def check_radianmore(self):
        rad1 = radian(90,1,60)
        assert_almost_equal(rad1,pi/2+0.0005816135199345904,5)

class test_reshape(unittest.TestCase):

    def check_reshape(self):
        a = (array([1,2,3]),array([4,5,6]))
        b = reshape(a,(2,3))
        assert_array_equal(b,array([[1, 2, 3],
                                    [4, 5, 6]]))
        c = reshape(a,(3,2))
        assert_array_equal(c,array([[1, 2],
                                    [3, 4],
                                    [5, 6]]))

class test_rgamma(unittest.TestCase):

    def check_rgamma(self):
        rgam = rgamma(8)
        rlgam = 1/gamma(8)
        assert_almost_equal(rgam,rlgam,8)

class test_riccati_jn(unittest.TestCase):

    def check_riccati_jn(self):
        jnrl = (sph_jn(1,.2)[0]*.2,sph_jn(1,.2)[0]+sph_jn(1,.2)[1]*.2)
        ricjn = riccati_jn(1,.2)
        assert_array_almost_equal(ricjn,jnrl,8)

class test_riccati_yn(unittest.TestCase):

    def check_riccati_yn(self):
        ynrl = (sph_yn(1,.2)[0]*.2,sph_yn(1,.2)[0]+sph_yn(1,.2)[1]*.2)
        ricyn = riccati_yn(1,.2)
        assert_array_almost_equal(ricyn,ynrl,8)

class test_round(unittest.TestCase):

    def check_round(self):
        rnd = map(int,(round(10.1),round(10.4),round(10.5),round(10.6)))
        rndrl = (10,10,11,11)
        assert_array_equal(rnd,rndrl)

class test_sh_legendre(unittest.TestCase):

    def check_sh_legendre(self):
        # P*_n(x) = P_n(2x-1)
        psub = poly1d([2,-1])
        Ps0 = sh_legendre(0)
        Ps1 = sh_legendre(1)        
        Ps2 = sh_legendre(2)
        Ps3 = sh_legendre(3)
        Ps4 = sh_legendre(4)
        Ps5 = sh_legendre(5)
        pse0 = legendre(0)(psub)
        pse1 = legendre(1)(psub)
        pse2 = legendre(2)(psub)
        pse3 = legendre(3)(psub)
        pse4 = legendre(4)(psub)
        pse5 = legendre(5)(psub)
        assert_array_almost_equal(Ps0.c,pse0.c,13)
        assert_array_almost_equal(Ps1.c,pse1.c,13)
        assert_array_almost_equal(Ps2.c,pse2.c,13)
        assert_array_almost_equal(Ps3.c,pse3.c,13)
        assert_array_almost_equal(Ps4.c,pse4.c,12)
        assert_array_almost_equal(Ps5.c,pse5.c,12)
        
class _test_sh_chebyt(unittest.TestCase):

    def check_sh_chebyt(self):
        # T*_n(x) = T_n(2x-1)
        psub = poly1d([2,-1])
        Ts0 = sh_chebyt(0)
        Ts1 = sh_chebyt(1)        
        Ts2 = sh_chebyt(2)
        Ts3 = sh_chebyt(3)
        Ts4 = sh_chebyt(4)
        Ts5 = sh_chebyt(5)
        tse0 = chebyt(0)(psub)
        tse1 = chebyt(1)(psub)
        tse2 = chebyt(2)(psub)
        tse3 = chebyt(3)(psub)
        tse4 = chebyt(4)(psub)
        tse5 = chebyt(5)(psub)
        assert_array_almost_equal(Ts0.c,tse0.c,13)
        assert_array_almost_equal(Ts1.c,tse1.c,13)
        assert_array_almost_equal(Ts2.c,tse2.c,13)
        assert_array_almost_equal(Ts3.c,tse3.c,13)
        assert_array_almost_equal(Ts4.c,tse4.c,12)
        assert_array_almost_equal(Ts5.c,tse5.c,12)
        

class test_sh_chebyu(unittest.TestCase):

    def check_sh_chebyu(self):
        # U*_n(x) = U_n(2x-1)
        psub = poly1d([2,-1])
        Us0 = sh_chebyu(0)
        Us1 = sh_chebyu(1)        
        Us2 = sh_chebyu(2)
        Us3 = sh_chebyu(3)
        Us4 = sh_chebyu(4)
        Us5 = sh_chebyu(5)
        use0 = chebyu(0)(psub)
        use1 = chebyu(1)(psub)
        use2 = chebyu(2)(psub)
        use3 = chebyu(3)(psub)
        use4 = chebyu(4)(psub)
        use5 = chebyu(5)(psub)
        assert_array_almost_equal(Us0.c,use0.c,13)
        assert_array_almost_equal(Us1.c,use1.c,13)
        assert_array_almost_equal(Us2.c,use2.c,13)
        assert_array_almost_equal(Us3.c,use3.c,13)
        assert_array_almost_equal(Us4.c,use4.c,12)
        assert_array_almost_equal(Us5.c,use5.c,11)

class test_sh_jacobi(unittest.TestCase):

    def check_sh_jacobi(self):
        # G^(p,q)_n(x) = n! gamma(n+p)/gamma(2*n+p) * P^(p-q,q-1)_n(2*x-1)
        conv = lambda n,p: _gam(n+1)*_gam(n+p)/_gam(2*n+p)
        psub = poly1d([2,-1])
        q = 4*rand()
        p = q-1 + 2*rand()
        print "shifted jacobi p,q = ", p, q
        G0 = sh_jacobi(0,p,q)
        G1 = sh_jacobi(1,p,q)        
        G2 = sh_jacobi(2,p,q)
        G3 = sh_jacobi(3,p,q)
        G4 = sh_jacobi(4,p,q)
        G5 = sh_jacobi(5,p,q)
        ge0 = jacobi(0,p-q,q-1)(psub) * conv(0,p)
        ge1 = jacobi(1,p-q,q-1)(psub) * conv(1,p)
        ge2 = jacobi(2,p-q,q-1)(psub) * conv(2,p)
        ge3 = jacobi(3,p-q,q-1)(psub) * conv(3,p)
        ge4 = jacobi(4,p-q,q-1)(psub) * conv(4,p)
        ge5 = jacobi(5,p-q,q-1)(psub) * conv(5,p)

        assert_array_almost_equal(G0.c,ge0.c,13)
        assert_array_almost_equal(G1.c,ge1.c,13)
        assert_array_almost_equal(G2.c,ge2.c,13)
        assert_array_almost_equal(G3.c,ge3.c,13)
        assert_array_almost_equal(G4.c,ge4.c,13)
        assert_array_almost_equal(G5.c,ge5.c,13)
        

class test_sinc(unittest.TestCase):

    def check_sinc(self):
        c = arange(-2,2,.1)
        y = sinc(c)
        yre = sin(pi*c)/(pi*c)
	yre[20] = 1.0
        assert_array_almost_equal(y, yre, 4)
    def check_0(self):
        x = 0.0
        assert_equal(sinc(x),1.0)

class test_sindg(unittest.TestCase):

    def check_sindg(self):
        sn = sindg(90)
        assert_equal(sn,1.0)

    def check_sindgmore(self):
        snm = sindg(30)
        snmrl = sin(pi/6.0)
        assert_almost_equal(snm,snmrl,8)
        snm1 = sindg(45)
        snmrl1 = sin(pi/4.0)
        assert_almost_equal(snm1,snmrl1,8)

class test_sph_harm(unittest.TestCase):

    def check_sph_harm(self):
        pass


class test_sph_in(unittest.TestCase):

    def check_sph_in(self):
        i1n = sph_in(1,.2)
        inp0 = (i1n[0][1])
        inp1 = (i1n[0][0] - 2.0/0.2 * i1n[0][1])
        assert_array_almost_equal(i1n[0],array([1.0066800127054699381,
                                                0.066933714568029540839]),12)
        assert_array_almost_equal(i1n[1],[inp0,inp1],12)

class test_sph_inkn(unittest.TestCase):

    def check_sph_inkn(self):
        spikn = r_[sph_in(1,.2)+sph_kn(1,.2)]
        inkn = r_[sph_inkn(1,.2)]
        assert_array_almost_equal(inkn,spikn,10)

class test_sph_jn(unittest.TestCase):

    def check_sph_jn(self):
        s1 = sph_jn(2,.2)
        s10 = -s1[0][1]
        s11 = s1[0][0]-2.0/0.2*s1[0][1]
        s12 = s1[0][1]-3.0/0.2*s1[0][2]
        assert_array_almost_equal(s1[0],[0.99334665397530607731,
                                      0.066400380670322230863,
                                      0.0026590560795273856680],12)
        assert_array_almost_equal(s1[1],[s10,s11,s12],12)

class test_sph_jnyn(unittest.TestCase):

    def check_sph_jnyn(self):
        jnyn = r_[sph_jn(1,.2) + sph_yn(1,.2)]  # tuple addition
        jnyn1 = r_[sph_jnyn(1,.2)]
        assert_array_almost_equal(jnyn1,jnyn,9)

class test_sph_kn(unittest.TestCase):

    def check_sph_kn(self):
        kn = sph_kn(2,.2)
        kn0 = -kn[0][1]
        kn1 = -kn[0][0]-2.0/0.2*kn[0][1]
        kn2 = -kn[0][1]-3.0/0.2*kn[0][2]
        assert_array_almost_equal(kn[0],[6.4302962978445670140,
                                         38.581777787067402086,
                                         585.15696310385559829],12)
        assert_array_almost_equal(kn[1],[kn0,kn1,kn2],9)

class test_sph_yn(unittest.TestCase):

    def check_sph_yn(self):
        sy1 = sph_yn(2,.2)[0][2]
        sy2 = sph_yn(0,.2)[0][0]
        sphpy = (sph_yn(1,.2)[0][0]-2*sph_yn(2,.2)[0][2])/3 #correct derivative value
        assert_almost_equal(sy1,-377.52483,5)#previous values in the system
        assert_almost_equal(sy2,-4.9003329,5)
        sy3 = sph_yn(1,.2)[1][1]
        assert_almost_equal(sy3,sphpy,4) #compare correct derivative val. (correct =-system val).

class test_take(unittest.TestCase):

    def check_take(self):
        a = array([0,1,2,3,4,5,6,7,8])
        tka = take(a,(0,4,5,8))
        assert_array_equal(tka,array([0,4,5,8]))

class test_tandg(unittest.TestCase):

    def check_tandg(self):
        tn = tandg(30)
        tnrl = tan(pi/6.0)
        assert_almost_equal(tn,tnrl,8)

    def check_tandgmore(self):
        tnm = tandg(45)
        tnmrl = tan(pi/4.0)
        assert_almost_equal(tnm,tnmrl,8)
        tnm1 = tandg(60)
        tnmrl1 = tan(pi/3.0)
        assert_almost_equal(tnm1,tnmrl1,8)

class test_y0(unittest.TestCase):

    def check_y0(self):
        oz = y0(.1)
        ozr = yn(0,.1)
        assert_almost_equal(oz,ozr,8)

class test_y1(unittest.TestCase):

    def check_y1(self):
        o1 = y1(.1)
        o1r = yn(1,.1)
        assert_almost_equal(o1,o1r,8)

class test_y0_zeros(unittest.TestCase):

    def check_y0_zeros(self):
        yo,ypo = y0_zeros(2)
        zo,zpo = y0_zeros(2,complex=1)
        all = r_[yo,zo]
        allval = r_[ypo,zpo]
        assert_array_almost_equal(abs(yv(0.0,all)),0.0,11)
        assert_array_almost_equal(abs(yv(1,all)-allval),0.0,11)
                                         

class test_y1_zeros(unittest.TestCase):

    def check_y1_zeros(self):
        y1 = y1_zeros(1)
        assert_array_almost_equal(y1,(array([2.19714]),array([0.52079])),5)
        
class test_y1p_zeros(unittest.TestCase):

    def check_y1p_zeros(self):
        y1p = y1p_zeros(1,complex=1)
        assert_array_almost_equal(y1p,(array([ 0.5768+0.904j]), array([-0.7635+0.5892j])),3)

class test_yn_zeros(unittest.TestCase):

    def check_yn_zeros(self):
        an = yn_zeros(4,2)
        assert_array_almost_equal(an,array([ 5.64515,  9.36162]),5)

class test_ynp_zeros(unittest.TestCase):

    def check_ynp_zeros(self):
        ao = ynp_zeros(0,2)
        assert_array_almost_equal(ao,array([ 2.19714133, 5.42968104]),6)

class test_yn(unittest.TestCase):

    def check_yn(self):
        yn2n = yn(1,.2)
        assert_almost_equal(yn2n,-3.3238249881118471,8)

class test_yv(unittest.TestCase):

    def check_yv(self):
        yv2 = yv(1,.2)
        assert_almost_equal(yv2,-3.3238249881118471,8)

class test_yve(unittest.TestCase):

    def check_yve(self):
        yve2 = yve(1,.2)
        assert_almost_equal(yve2,-3.3238249881118471,8)
        yve2r = yv(1,.2+1j)*exp(-1)
        yve22 = yve(1,.2+1j)
        assert_almost_equal(yve22,yve2r,8)

class test_yvp(unittest.TestCase):

    def check_yvp(self):
        yvpr = (yv(1,.2) - yv(3,.2))/2.0
        yvp1 = yvp(2,.2)
        assert_array_almost_equal(yvp1,yvpr,6)

class test_zeros(unittest.TestCase):

    def check_zeros(self):
        b = zeros((1,11))
        assert_array_equal(b,array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]))
        c = zeros((2,2))
        assert_array_equal(c,array([[0, 0],
                                    [0, 0]]))

if __name__ == "__main__":
    ScipyTest('special.basic').run()
