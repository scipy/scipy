#this program corresponds to special.py

### Means test is not done yet
#E   Means test is giving error (E)
#F   Means test is failing (F)
#EF  Means test is giving error and Failing
#!   Means test is segfaulting
#8   Means test runs forever

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
#8   test_sh_chebyu
#8   test_sh_jacobi
#8   test_sh_legendre

from numpy import dot, array

from scipy.testing import *

from scipy.special import *
import scipy.special._cephes as cephes


class TestCephes(TestCase):
    def test_airy(self):
        cephes.airy(0)
    def test_airye(self):
        cephes.airye(0)
    def test_bdtr(self):
        assert_equal(cephes.bdtr(1,1,0.5),1.0)
    def test_bdtri(self):
        assert_equal(cephes.bdtri(1,3,0.5),0.5)
    def test_bdtrc(self):
        assert_equal(cephes.bdtrc(1,3,0.5),0.5)
    def test_bdtrin(self):
        assert_equal(cephes.bdtrin(1,0,1),5.0)
    def test_bdtrik(self):
        cephes.bdtrik(1,3,0.5)

    def test_bei(self):
        assert_equal(cephes.bei(0),0.0)
    def test_beip(self):
        assert_equal(cephes.beip(0),0.0)
    def test_ber(self):
        assert_equal(cephes.ber(0),1.0)
    def test_berp(self):
        assert_equal(cephes.berp(0),0.0)

    def test_besselpoly(self):
        assert_equal(cephes.besselpoly(0,0,0),1.0)

    def test_beta(self):
        assert_equal(cephes.beta(1,1),1.0)
    def test_betainc(self):
        assert_equal(cephes.betainc(1,1,1),1.0)
    def test_betaln(self):
        assert_equal(cephes.betaln(1,1),0.0)
    def test_betaincinv(self):
        assert_equal(cephes.betaincinv(1,1,1),1.0)

    def test_btdtr(self):
        assert_equal(cephes.btdtr(1,1,1),1.0)
    def test_btdtri(self):
        assert_equal(cephes.btdtri(1,1,1),1.0)
    def test_btdtria(self):
        assert_equal(cephes.btdtria(1,1,1),5.0)
    def test_btdtrib(self):
        assert_equal(cephes.btdtrib(1,1,1),5.0)

    def test_cbrt(self):
        assert_approx_equal(cephes.cbrt(1),1.0)

    def test_chdtr(self):
        assert_equal(cephes.chdtr(1,0),0.0)
    def test_chdtrc(self):
        assert_equal(cephes.chdtrc(1,0),1.0)
    def test_chdtri(self):
        assert_equal(cephes.chdtri(1,1),0.0)
    def test_chdtriv(self):
        assert_equal(cephes.chdtriv(0,0),5.0)

    def test_chndtr(self):
        assert_equal(cephes.chndtr(0,1,0),0.0)
    def test_chndtridf(self):
        assert_equal(cephes.chndtridf(0,0,1),5.0)
    def test_chndtrinc(self):
        assert_equal(cephes.chndtrinc(0,1,0),5.0)
    def test_chndtrix(self):
        assert_equal(cephes.chndtrix(0,1,0),0.0)

    def test_cosdg(self):
        assert_equal(cephes.cosdg(0),1.0)
    def test_cosm1(self):
        assert_equal(cephes.cosm1(0),0.0)
    def test_cotdg(self):
        assert_almost_equal(cephes.cotdg(45),1.0)

    def test_dawsn(self):
        assert_equal(cephes.dawsn(0),0.0)

    def test_ellipe(self):
        assert_equal(cephes.ellipe(1),1.0)
    def test_ellipeinc(self):
        assert_equal(cephes.ellipeinc(0,1),0.0)
    def test_ellipj(self):
        cephes.ellipj(0,1)
    def test_ellipk(self):
        cephes.ellipk(0)#==pi/2
    def test_ellipkinc(self):
        assert_equal(cephes.ellipkinc(0,0),0.0)

    def test_erf(self):
        assert_equal(cephes.erf(0),0.0)
    def test_erfc(self):
        assert_equal(cephes.erfc(0),1.0)

    def test_exp1(self):
        cephes.exp1(1)
    def test_expi(self):
        cephes.expi(1)
    def test_expn(self):
        cephes.expn(1,1)

    def test_exp10(self):
        assert_approx_equal(cephes.exp10(2),100.0)
    def test_exp2(self):
        assert_equal(cephes.exp2(2),4.0)
    def test_expm1(self):
        assert_equal(cephes.expm1(0),0.0)

    def test_fdtr(self):
        assert_equal(cephes.fdtr(1,1,0),0.0)
    def test_fdtrc(self):
        assert_equal(cephes.fdtrc(1,1,0),1.0)
    def test_fdtri(self):
        cephes.fdtri(1,1,0.5)
    def test_fdtridfd(self):
        assert_equal(cephes.fdtridfd(1,0,0),5.0)

    def test_fresnel(self):
        assert_equal(cephes.fresnel(0),(0.0,0.0))

    def test_gamma(self):
        assert_equal(cephes.gamma(5),24.0)
    def test_gammainc(self):
        assert_equal(cephes.gammainc(5,0),0.0)
    def test_gammaincc(self):
        assert_equal(cephes.gammaincc(5,0),1.0)
    def test_gammainccinv(self):
        assert_equal(cephes.gammainccinv(5,1),0.0)
    def test_gammaln(self):
        cephes.gammaln(10)

    def test_gdtr(self):
        assert_equal(cephes.gdtr(1,1,0),0.0)
    def test_gdtrc(self):
        assert_equal(cephes.gdtrc(1,1,0),1.0)
    def test_gdtria(self):
        assert_equal(cephes.gdtria(0,1,1),0.0)
    def test_gdtrib(self):
        cephes.gdtrib(1,0,1)
        #assert_equal(cephes.gdtrib(1,0,1),5.0)
    def test_gdtrix(self):
        cephes.gdtrix(1,1,.1)

    def test_hankel1(self):
        cephes.hankel1(1,1)
    def test_hankel1e(self):
        cephes.hankel1e(1,1)
    def test_hankel2(self):
        cephes.hankel2(1,1)
    def test_hankel2e(self):
        cephes.hankel2e(1,1)

    def test_hyp1f1(self):
        assert_approx_equal(cephes.hyp1f1(1,1,1), exp(1.0))
        assert_approx_equal(cephes.hyp1f1(3,4,-6), 0.026056422099537251095)
        cephes.hyp1f1(1,1,1)
    def test_hyp1f2(self):
        cephes.hyp1f2(1,1,1,1)
    def test_hyp2f0(self):
        cephes.hyp2f0(1,1,1,1)
    def test_hyp2f1(self):
        assert_equal(cephes.hyp2f1(1,1,1,0),1.0)
    def test_hyp3f0(self):
        assert_equal(cephes.hyp3f0(1,1,1,0),(1.0,0.0))
    def test_hyperu(self):
        assert_equal(cephes.hyperu(0,1,1),1.0)

    def test_i0(self):
        assert_equal(cephes.i0(0),1.0)
    def test_i0e(self):
        assert_equal(cephes.i0e(0),1.0)
    def test_i1(self):
        assert_equal(cephes.i1(0),0.0)
    def test_i1e(self):
        assert_equal(cephes.i1e(0),0.0)

    def test_it2i0k0(self):
        cephes.it2i0k0(1)
    def test_it2j0y0(self):
        cephes.it2j0y0(1)
    def test_it2struve0(self):
        cephes.it2struve0(1)
    def test_itairy(self):
        cephes.itairy(1)
    def test_iti0k0(self):
        assert_equal(cephes.iti0k0(0),(0.0,0.0))
    def test_itj0y0(self):
        assert_equal(cephes.itj0y0(0),(0.0,0.0))
    def test_itmodstruve0(self):
        assert_equal(cephes.itmodstruve0(0),0.0)
    def test_itstruve0(self):
        assert_equal(cephes.itstruve0(0),0.0)
    def test_iv(self):
        assert_equal(cephes.iv(1,0),0.0)
    def _check_ive(self):
        assert_equal(cephes.ive(1,0),0.0)

    def test_j0(self):
        assert_equal(cephes.j0(0),1.0)
    def test_j1(self):
        assert_equal(cephes.j1(0),0.0)
    def test_jn(self):
        assert_equal(cephes.jn(0,0),1.0)
    def test_jv(self):
        assert_equal(cephes.jv(0,0),1.0)
    def _check_jve(self):
        assert_equal(cephes.jve(0,0),1.0)

    def test_k0(self):
        cephes.k0(2)
    def test_k0e(self):
        cephes.k0e(2)
    def test_k1(self):
        cephes.k1(2)
    def test_k1e(self):
        cephes.k1e(2)
    def test_kei(self):
        cephes.kei(2)
    def test_keip(self):
        assert_equal(cephes.keip(0),0.0)
    def test_ker(self):
        cephes.ker(2)
    def test_kerp(self):
        cephes.kerp(2)
    def _check_kelvin(self):
        cephes.kelvin(2)
    def test_kn(self):
        cephes.kn(1,1)

    def test_kolmogi(self):
        assert_equal(cephes.kolmogi(1),0.0)
    def test_kolmogorov(self):
        assert_equal(cephes.kolmogorov(0),1.0)

    def _check_kv(self):
        cephes.kv(1,1)
    def _check_kve(self):
        cephes.kve(1,1)
    def test_log1p(self):
        assert_equal(cephes.log1p(0),0.0)
    def test_lpmv(self):
        assert_equal(cephes.lpmv(0,0,1),1.0)

    def test_mathieu_a(self):
        assert_equal(cephes.mathieu_a(1,0),1.0)
    def test_mathieu_b(self):
        assert_equal(cephes.mathieu_b(1,0),1.0)
    def test_mathieu_cem(self):
        assert_equal(cephes.mathieu_cem(1,0,0),(1.0,0.0))
    def test_mathieu_modcem1(self):
        assert_equal(cephes.mathieu_modcem1(1,0,0),(0.0,0.0))
    def test_mathieu_modcem2(self):
        cephes.mathieu_modcem2(1,1,1)
    def test_mathieu_sem(self):
        assert_equal(cephes.mathieu_sem(1,0,0),(0.0,1.0))
    def test_mathieu_modsem1(self):
        assert_equal(cephes.mathieu_modsem1(1,0,0),(0.0,0.0))
    def test_mathieu_modsem2(self):
        cephes.mathieu_modsem2(1,1,1)

    def test_modfresnelm(self):
        cephes.modfresnelm(0)
    def test_modfresnelp(self):
        cephes.modfresnelp(0)
    def _check_modstruve(self):
        assert_equal(cephes.modstruve(1,0),0.0)

    def test_nbdtr(self):
        assert_equal(cephes.nbdtr(1,1,1),1.0)
    def test_nbdtrc(self):
        assert_equal(cephes.nbdtrc(1,1,1),0.0)
    def test_nbdtri(self):
        assert_equal(cephes.nbdtri(1,1,1),1.0)
    def __check_nbdtrik(self):
        cephes.nbdtrik(1,.4,.5)
    def test_nbdtrin(self):
        assert_equal(cephes.nbdtrin(1,0,0),5.0)

    def test_ncfdtr(self):
        assert_equal(cephes.ncfdtr(1,1,1,0),0.0)
    def test_ncfdtri(self):
        assert_equal(cephes.ncfdtri(1,1,1,0),0.0)
    def test_ncfdtridfd(self):
        cephes.ncfdtridfd(1,0.5,0,1)
    def __check_ncfdtridfn(self):
        cephes.ncfdtridfn(1,0.5,0,1)
    def __check_ncfdtrinc(self):
        cephes.ncfdtrinc(1,0.5,0,1)

    def test_nctdtr(self):
        assert_equal(cephes.nctdtr(1,0,0),0.5)
    def __check_nctdtridf(self):
        cephes.nctdtridf(1,0.5,0)
    def test_nctdtrinc(self):
        cephes.nctdtrinc(1,0,0)
    def test_nctdtrit(self):
        cephes.nctdtrit(.1,0.2,.5)

    def test_ndtr(self):
        assert_equal(cephes.ndtr(0),0.5)
    def test_ndtri(self):
        assert_equal(cephes.ndtri(0.5),0.0)
    def test_nrdtrimn(self):
        assert_approx_equal(cephes.nrdtrimn(0.5,1,1),1.0)
    def test_nrdtrisd(self):
        assert_equal(cephes.nrdtrisd(0.5,0.5,0.5),0.0)

    def test_obl_ang1(self):
        cephes.obl_ang1(1,1,1,0)
    def test_obl_ang1_cv(self):
        result = cephes.obl_ang1_cv(1,1,1,1,0)
        assert_almost_equal(result[0],1.0)
        assert_almost_equal(result[1],0.0)

    def _check_obl_cv(self):
        assert_equal(cephes.obl_cv(1,1,0),2.0)
    def test_obl_rad1(self):
        cephes.obl_rad1(1,1,1,0)
    def test_obl_rad1_cv(self):
        cephes.obl_rad1_cv(1,1,1,1,0)
    def test_obl_rad2(self):
        cephes.obl_rad2(1,1,1,0)
    def test_obl_rad2_cv(self):
        cephes.obl_rad2_cv(1,1,1,1,0)

    def test_pbdv(self):
        assert_equal(cephes.pbdv(1,0),(0.0,0.0))
    def test_pbvv(self):
        cephes.pbvv(1,0)
    def test_pbwa(self):
        cephes.pbwa(1,0)
    def test_pdtr(self):
        cephes.pdtr(0,1)
    def test_pdtrc(self):
        cephes.pdtrc(0,1)
    def test_pdtri(self):
        cephes.pdtri(0.5,0.5)
    def test_pdtrik(self):
        cephes.pdtrik(0.5,1)

    def test_pro_ang1(self):
        cephes.pro_ang1(1,1,1,0)
    def test_pro_ang1_cv(self):
        assert_array_almost_equal(cephes.pro_ang1_cv(1,1,1,1,0),
                                  array((1.0,0.0)))
    def _check_pro_cv(self):
        assert_equal(cephes.pro_cv(1,1,0),2.0)
    def test_pro_rad1(self):
        cephes.pro_rad1(1,1,1,0.1)
    def test_pro_rad1_cv(self):
        cephes.pro_rad1_cv(1,1,1,1,0)
    def test_pro_rad2(self):
        cephes.pro_rad2(1,1,1,0)
    def test_pro_rad2_cv(self):
        cephes.pro_rad2_cv(1,1,1,1,0)

    def test_psi(self):
        cephes.psi(1)

    def test_radian(self):
        assert_equal(cephes.radian(0,0,0),0)
    def test_rgamma(self):
        assert_equal(cephes.rgamma(1),1.0)
    def test_round(self):
        assert_equal(cephes.round(3.4),3.0)
        assert_equal(cephes.round(-3.4),-3.0)
        assert_equal(cephes.round(3.6),4.0)
        assert_equal(cephes.round(-3.6),-4.0)
        assert_equal(cephes.round(3.5),4.0)
        assert_equal(cephes.round(-3.5),-4.0)

    def test_shichi(self):
        cephes.shichi(1)
    def test_sici(self):
        cephes.sici(1)
    def test_sindg(self):
        assert_equal(cephes.sindg(90),1.0)
    def test_smirnov(self):
        assert_equal(cephes.smirnov(1,.1),0.9)
    def test_smirnovi(self):
        assert_almost_equal(cephes.smirnov(1,cephes.smirnovi(1,0.4)),0.4)
        assert_almost_equal(cephes.smirnov(1,cephes.smirnovi(1,0.6)),0.6)

    def test_spence(self):
        assert_equal(cephes.spence(1),0.0)
    def test_stdtr(self):
        assert_equal(cephes.stdtr(1,0),0.5)
    def test_stdtridf(self):
        cephes.stdtridf(0.7,1)
    def test_stdtrit(self):
        cephes.stdtrit(1,0.7)
    def test_struve(self):
        assert_equal(cephes.struve(0,0),0.0)

    def test_tandg(self):
        assert_equal(cephes.tandg(45),1.0)
    def test_tklmbda(self):
        assert_almost_equal(cephes.tklmbda(1,1),1.0)

    def test_y0(self):
        cephes.y0(1)
    def test_y1(self):
        cephes.y1(1)
    def test_yn(self):
        cephes.yn(1,1)
    def test_yv(self):
        cephes.yv(1,1)
    def _check_yve(self):
        cephes.yve(1,1)

    def test_zeta(self):
        cephes.zeta(2,2)
    def test_zetac(self):
        assert_equal(cephes.zetac(0),-1.5)
    def test_wofz(self):
        cephes.wofz(0)

class TestAiry(TestCase):
    def test_airy(self):
        #This tests the airy function to ensure 8 place accuracy in computation

        x = airy(.99)
        assert_array_almost_equal(x,array([0.13689066,-0.16050153,1.19815925,0.92046818]),8)
        x = airy(.41)
        assert_array_almost_equal(x,array([0.25238916,-.23480512,0.80686202,0.51053919]),8)
        x = airy(-.36)
        assert_array_almost_equal(x,array([0.44508477,-0.23186773,0.44939534,0.48105354]),8)

    def test_airye(self):
        a = airye(0.01)
        b = airy(0.01)
        b1 = [None]*4
        for n in range(2):
            b1[n] = b[n]*exp(2.0/3.0*0.01*sqrt(0.01))
        for n in range(2,4):
            b1[n] = b[n]*exp(-abs(real(2.0/3.0*0.01*sqrt(0.01))))
        assert_array_almost_equal(a,b1,6)

    def test_bi_zeros(self):
        bi = bi_zeros(2)
        bia = (array([-1.17371322, -3.2710930]),
        array([-2.29443968, -4.07315509]),
        array([-0.45494438,  0.39652284]),
        array([ 0.60195789 , -0.76031014]))
        assert_array_almost_equal(bi,bia,4)

    def test_ai_zeros(self):
        ai = ai_zeros(1)
        assert_array_almost_equal(ai,(array([-2.33810741]),
                                     array([-1.01879297]),
                                     array([ 0.5357]),
                                     array([ 0.7012])),4)

class TestAssocLaguerre(TestCase):
    def test_assoc_laguerre(self):
        a1 = genlaguerre(11,1)
        a2 = assoc_laguerre(.2,11,1)
        assert_array_almost_equal(a2,a1(.2),8)
        a2 = assoc_laguerre(1,11,1)
        assert_array_almost_equal(a2,a1(1),8)

class TestBesselpoly(TestCase):
    def test_besselpoly(self):
        pass

class TestKelvin(TestCase):
    def test_bei(self):
        mbei = bei(2)
        assert_almost_equal(mbei, 0.9722916273066613,5)#this may not be exact

    def test_beip(self):
        mbeip = beip(2)
        assert_almost_equal(mbeip,0.91701361338403631,5)#this may not be exact

    def test_ber(self):
        mber = ber(2)
        assert_almost_equal(mber,0.75173418271380821,5)#this may not be exact

    def test_berp(self):
        mberp = berp(2)
        assert_almost_equal(mberp,-0.49306712470943909,5)#this may not be exact

    def test_bei_zeros(self):
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


    def test_beip_zeros(self):
        bip = beip_zeros(5)
        assert_array_almost_equal(bip,array([  3.772673304934953,
                                               8.280987849760042,
                                               12.742147523633703,
                                               17.193431752512542,
                                               21.641143941167325]),4)

    def test_ber_zeros(self):
        ber = ber_zeros(5)
        assert_array_almost_equal(ber,array([2.84892,
                                             7.23883,
                                             11.67396,
                                             16.11356,
                                             20.55463]),4)

    def test_berp_zeros(self):
        brp = berp_zeros(5)
        assert_array_almost_equal(brp,array([6.03871,
                                             10.51364,
                                             14.96844,
                                             19.41758,
                                             23.86430]),4)

    def test_kelvin(self):
        mkelv = kelvin(2)
        assert_array_almost_equal(mkelv,(ber(2)+bei(2)*1j,
                                         ker(2)+kei(2)*1j,
                                         berp(2)+beip(2)*1j,
                                         kerp(2)+keip(2)*1j),8)

    def test_kei(self):
        mkei = kei(2)
        assert_almost_equal(mkei,-0.20240006776470432,5)

    def test_keip(self):
        mkeip = keip(2)
        assert_almost_equal(mkeip,0.21980790991960536,5)

    def test_ker(self):
        mker = ker(2)
        assert_almost_equal(mker,-0.041664513991509472,5)

    def test_kerp(self):
        mkerp = kerp(2)
        assert_almost_equal(mkerp,-0.10660096588105264,5)

    def test_kei_zeros(self):
        kei = kei_zeros(5)
        assert_array_almost_equal(kei,array([  3.91467,
                                              8.34422,
                                              12.78256,
                                              17.22314,
                                              21.66464]),4)

    def test_keip_zeros(self):
        keip = keip_zeros(5)
        assert_array_almost_equal(keip,array([  4.93181,
                                                9.40405,
                                                13.85827,
                                                18.30717,
                                                22.75379]),4)



    # numbers come from 9.9 of A&S pg. 381
    def test_kelvin_zeros(self):
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
        assert_array_almost_equal(kerpz,array([ 2.66584,
                                                7.17212,
                                                11.63218,
                                                16.08312,
                                                20.53068]),4)
        assert_array_almost_equal(keipz,array([ 4.93181,
                                                9.40405,
                                                13.85827,
                                                18.30717,
                                                22.75379]),4)

    def test_ker_zeros(self):
        ker = ker_zeros(5)
        assert_array_almost_equal(ker,array([  1.71854,
                                               6.12728,
                                               10.56294,
                                               15.00269,
                                               19.44381]),4)

    def test_kerp_zeros(self):
        kerp = kerp_zeros(5)
        assert_array_almost_equal(kerp,array([  2.66584,
                                                7.17212,
                                                11.63218,
                                                16.08312,
                                                20.53068]),4)

class TestBernoulli(TestCase):
    def test_bernoulli(self):
        brn = bernoulli(5)
        assert_array_almost_equal(brn,array([1.0000,
                                             -0.5000,
                                             0.1667,
                                             0.0000,
                                             -0.0333,
                                             0.0000]),4)

class TestBeta(TestCase):
    def test_beta(self):
        bet = beta(2,4)
        betg = (gamma(2)*gamma(4))/gamma(6)
        assert_almost_equal(bet,betg,8)

    def test_betaln(self):
        betln = betaln(2,4)
        bet = log(abs(beta(2,4)))
        assert_almost_equal(betln,bet,8)

    def test_betainc(self):
        btinc = betainc(1,1,.2)
        assert_almost_equal(btinc,0.2,8)

    def test_betaincinv(self):
        y = betaincinv(2,4,.5)
        comp = betainc(2,4,y)
        assert_almost_equal(comp,.5,5)

class TestCheby(TestCase):
    def test_chebyc(self):
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

    def test_chebys(self):
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

    def test_chebyt(self):
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

    def test_chebyu(self):
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

class TestTrigonometric(TestCase):
    def test_cbrt(self):
        cb = cbrt(27)
        cbrl = 27**(1.0/3.0)
        assert_approx_equal(cb,cbrl)

    def test_cbrtmore(self):
        cb1 = cbrt(27.9)
        cbrl1 = 27.9**(1.0/3.0)
        assert_almost_equal(cb1,cbrl1,8)

    def test_cosdg(self):
        cdg = cosdg(90)
        cdgrl = cos(pi/2.0)
        assert_almost_equal(cdg,cdgrl,8)

    def test_cosdgmore(self):
        cdgm = cosdg(30)
        cdgmrl = cos(pi/6.0)
        assert_almost_equal(cdgm,cdgmrl,8)

    def test_cosm1(self):
        cs = (cosm1(0),cosm1(.3),cosm1(pi/10))
        csrl = (cos(0)-1,cos(.3)-1,cos(pi/10)-1)
        assert_array_almost_equal(cs,csrl,8)

    def test_cotdg(self):
        ct = cotdg(30)
        ctrl = tan(pi/6.0)**(-1)
        assert_almost_equal(ct,ctrl,8)

    def test_cotdgmore(self):
        ct1 = cotdg(45)
        ctrl1 = tan(pi/4.0)**(-1)
        assert_almost_equal(ct1,ctrl1,8)

    def test_specialpoints(self):
        assert_almost_equal(cotdg(45), 1.0, 14)
        assert_almost_equal(cotdg(-45), -1.0, 14)
        assert_almost_equal(cotdg(90), 0.0, 14)
        assert_almost_equal(cotdg(-90), 0.0, 14)
        assert_almost_equal(cotdg(135), -1.0, 14)
        assert_almost_equal(cotdg(-135), 1.0, 14)
        assert_almost_equal(cotdg(225), 1.0, 14)
        assert_almost_equal(cotdg(-225), -1.0, 14)
        assert_almost_equal(cotdg(270), 0.0, 14)
        assert_almost_equal(cotdg(-270), 0.0, 14)
        assert_almost_equal(cotdg(315), -1.0, 14)
        assert_almost_equal(cotdg(-315), 1.0, 14)
        assert_almost_equal(cotdg(765), 1.0, 14)

    def test_sinc(self):
        c = arange(-2,2,.1)
        y = sinc(c)
        yre = sin(pi*c)/(pi*c)
        yre[20] = 1.0
        assert_array_almost_equal(y, yre, 4)
    def test_0(self):
        x = 0.0
        assert_equal(sinc(x),1.0)

    def test_sindg(self):
        sn = sindg(90)
        assert_equal(sn,1.0)

    def test_sindgmore(self):
        snm = sindg(30)
        snmrl = sin(pi/6.0)
        assert_almost_equal(snm,snmrl,8)
        snm1 = sindg(45)
        snmrl1 = sin(pi/4.0)
        assert_almost_equal(snm1,snmrl1,8)

class TestTandg(TestCase):

    def test_tandg(self):
        tn = tandg(30)
        tnrl = tan(pi/6.0)
        assert_almost_equal(tn,tnrl,8)

    def test_tandgmore(self):
        tnm = tandg(45)
        tnmrl = tan(pi/4.0)
        assert_almost_equal(tnm,tnmrl,8)
        tnm1 = tandg(60)
        tnmrl1 = tan(pi/3.0)
        assert_almost_equal(tnm1,tnmrl1,8)

    def test_specialpoints(self):
        assert_almost_equal(tandg(0), 0.0, 14)
        assert_almost_equal(tandg(45), 1.0, 14)
        assert_almost_equal(tandg(-45), -1.0, 14)
        assert_almost_equal(tandg(135), -1.0, 14)
        assert_almost_equal(tandg(-135), 1.0, 14)
        assert_almost_equal(tandg(180), 0.0, 14)
        assert_almost_equal(tandg(-180), 0.0, 14)
        assert_almost_equal(tandg(225), 1.0, 14)
        assert_almost_equal(tandg(-225), -1.0, 14)
        assert_almost_equal(tandg(315), -1.0, 14)
        assert_almost_equal(tandg(-315), 1.0, 14)

class TestEllip(TestCase):
    def test_ellipj(self):
        el = ellipj(0.2,0)
        rel = [sin(0.2),cos(0.2),1.0,0.20]
        assert_array_almost_equal(el,rel,13)

    def test_ellipk(self):
        elk = ellipk(.2)
        assert_almost_equal(elk,1.659623598610528,11)

    def test_ellipkinc(self):
        elkinc = ellipkinc(pi/2,.2)
        elk = ellipk(0.2)
        assert_almost_equal(elkinc,elk,15)
        alpha = 20*pi/180
        phi = 45*pi/180
        m = sin(alpha)**2
        elkinc = ellipkinc(phi,m)
        assert_almost_equal(elkinc,0.79398143,8)
        # From pg. 614 of A & S

    def test_ellipe(self):
        ele = ellipe(.2)
        assert_almost_equal(ele,1.4890350580958529,8)

    def test_ellipeinc(self):
        eleinc = ellipeinc(pi/2,.2)
        ele = ellipe(0.2)
        assert_almost_equal(eleinc,ele,14)
        # pg 617 of A & S
        alpha, phi = 52*pi/180,35*pi/180
        m = sin(alpha)**2
        eleinc = ellipeinc(phi,m)
        assert_almost_equal(eleinc, 0.58823065, 8)


class TestErf(TestCase):

    def test_erf(self):
        er = erf(.25)
        assert_almost_equal(er,0.2763263902,8)

    def test_erf_zeros(self):
        erz = erf_zeros(5)
        erzr= array([1.45061616+1.88094300j,
                     2.24465928+2.61657514j,
                     2.83974105+3.17562810j,
                     3.33546074+3.64617438j,
                     3.76900557+4.06069723j])
        assert_array_almost_equal(erz,erzr,4)

    def test_erfcinv(self):
        i = erfcinv(1)
        assert_equal(i,0)

    def test_erfinv(self):
        i = erfinv(0)
        assert_equal(i,0)

    def test_errprint(self):
        a = errprint()
        b = 1-a #a is the state 1-a inverts state
        c = errprint(b) #returns last state 'a'
        assert_equal(a,c)
        d = errprint(a) #returns to original state
        assert_equal(d,b) #makes sure state was returned
        #assert_equal(d,1-a)

class TestEuler(TestCase):
    def test_euler(self):
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

class TestExp(TestCase):
    def test_exp2(self):
        ex = exp2(2)
        exrl = 2**2
        assert_equal(ex,exrl)

    def test_exp2more(self):
        exm = exp2(2.5)
        exmrl = 2**(2.5)
        assert_almost_equal(exm,exmrl,8)

    def test_exp10(self):
        ex = exp10(2)
        exrl = 10**2
        assert_approx_equal(ex,exrl)

    def test_exp10more(self):
        exm = exp10(2.5)
        exmrl = 10**(2.5)
        assert_almost_equal(exm,exmrl,8)

    def test_expm1(self):
        ex = (expm1(2),expm1(3),expm1(4))
        exrl = (exp(2)-1,exp(3)-1,exp(4)-1)
        assert_array_almost_equal(ex,exrl,8)

    def test_expm1more(self):
        ex1 = (expm1(2),expm1(2.1),expm1(2.2))
        exrl1 = (exp(2)-1,exp(2.1)-1,exp(2.2)-1)
        assert_array_almost_equal(ex1,exrl1,8)

class TestFresnel(TestCase):
    def test_fresnel(self):
        frs = array(fresnel(.5))
        assert_array_almost_equal(frs,array([0.064732432859999287, 0.49234422587144644]),8)

    # values from pg 329  Table 7.11 of A & S
    #  slightly corrected in 4th decimal place
    def test_fresnel_zeros(self):
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

    def test_fresnelc_zeros(self):
        szo, czo = fresnel_zeros(6)
        frc = fresnelc_zeros(6)
        assert_array_almost_equal(frc,czo,12)

    def test_fresnels_zeros(self):
        szo, czo = fresnel_zeros(5)
        frs = fresnels_zeros(5)
        assert_array_almost_equal(frs,szo,12)


class TestGamma(TestCase):
    def test_gamma(self):
        gam = gamma(5)
        assert_equal(gam,24.0)

    def test_gammaln(self):
        gamln = gammaln(3)
        lngam = log(gamma(3))
        assert_almost_equal(gamln,lngam,8)

    def test_gammainc(self):
        gama = gammainc(.5,.5)
        assert_almost_equal(gama,.7,1)

    def test_gammaincc(self):
        gicc = gammaincc(.5,.5)
        greal = 1 - gammainc(.5,.5)
        assert_almost_equal(gicc,greal,8)

    def test_gammainccinv(self):
        gccinv = gammainccinv(.5,.5)
        gcinv = gammaincinv(.5,.5)
        assert_almost_equal(gccinv,gcinv,8)

    def test_gammaincinv(self):
        y = gammaincinv(.4,.4)
        x = gammainc(.4,y)
        assert_almost_equal(x,0.4,1)
        y = gammainc(10, 0.05)
        x = gammaincinv(10, 2.5715803516000736e-20)
        assert_almost_equal(0.05, x, decimal=10)
        assert_almost_equal(y, 2.5715803516000736e-20, decimal=10)
        x = gammaincinv(50, 8.20754777388471303050299243573393e-18)
        assert_almost_equal(11.0, x, decimal=10)

    def test_rgamma(self):
        rgam = rgamma(8)
        rlgam = 1/gamma(8)
        assert_almost_equal(rgam,rlgam,8)

class TestHankel(TestCase):
    def test_negv(self):
        assert_almost_equal(hankel1(-3,2), -hankel1(3,2), 14)

    def test_hankel1(self):
        hank1 = hankel1(1,.1)
        hankrl = (jv(1,.1)+yv(1,.1)*1j)
        assert_almost_equal(hank1,hankrl,8)

    def test_negv(self):
        assert_almost_equal(hankel1e(-3,2), -hankel1e(3,2), 14)

    def test_hankel1e(self):
        hank1e = hankel1e(1,.1)
        hankrle = hankel1(1,.1)*exp(-.1j)
        assert_almost_equal(hank1e,hankrle,8)

    def test_negv(self):
        assert_almost_equal(hankel2(-3,2), -hankel2(3,2), 14)

    def test_hankel2(self):
        hank2 = hankel2(1,.1)
        hankrl2 = (jv(1,.1)-yv(1,.1)*1j)
        assert_almost_equal(hank2,hankrl2,8)

    def test_negv(self):
        assert_almost_equal(hankel2e(-3,2), -hankel2e(3,2), 14)

    def test_hankl2e(self):
        hank2e = hankel2e(1,.1)
        hankrl2e = hankel2e(1,.1)
        assert_almost_equal(hank2e,hankrl2e,8)

class TestHermite(TestCase):
    def test_hermite(self):
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

    def test_hermitenorm(self):
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

class TestGegenbauer(TestCase):

    def test_gegenbauer(self):
        a = 5*rand()-0.5
        if any(a==0): a = -0.2
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


class TestHyper(TestCase):
    def test_h1vp(self):
        h1 = h1vp(1,.1)
        h1real = (jvp(1,.1)+yvp(1,.1)*1j)
        assert_almost_equal(h1,h1real,8)

    def test_h2vp(self):
        h2 = h2vp(1,.1)
        h2real = (jvp(1,.1)-yvp(1,.1)*1j)
        assert_almost_equal(h2,h2real,8)

    def test_hyp0f1(self):
        pass

    def test_hyp1f1(self):
        hyp1 = hyp1f1(.1,.1,.3)
        assert_almost_equal(hyp1, 1.3498588075760032,7)

        # test contributed by Moritz Deger (2008-05-29)
        # http://projects.scipy.org/scipy/scipy/ticket/659
        
        # reference data obtained from mathematica [ a, b, x, m(a,b,x)]:
        # produced with test_hyp1f1.nb
        ref_data = array([[ -8.38132975e+00,  -1.28436461e+01,  -2.91081397e+01,          1.04178330e+04],
                          [  2.91076882e+00,  -6.35234333e+00,  -1.27083993e+01,          6.68132725e+00],
                          [ -1.42938258e+01,   1.80869131e-01,   1.90038728e+01,          1.01385897e+05],
                          [  5.84069088e+00,   1.33187908e+01,   2.91290106e+01,          1.59469411e+08],
                          [ -2.70433202e+01,  -1.16274873e+01,  -2.89582384e+01,          1.39900152e+24],
                          [  4.26344966e+00,  -2.32701773e+01,   1.91635759e+01,          6.13816915e+21],
                          [  1.20514340e+01,  -3.40260240e+00,   7.26832235e+00,          1.17696112e+13],
                          [  2.77372955e+01,  -1.99424687e+00,   3.61332246e+00,          3.07419615e+13],
                          [  1.50310939e+01,  -2.91198675e+01,  -1.53581080e+01,         -3.79166033e+02],
                          [  1.43995827e+01,   9.84311196e+00,   1.93204553e+01,          2.55836264e+10],
                          [ -4.08759686e+00,   1.34437025e+01,  -1.42072843e+01,          1.70778449e+01],
                          [  8.05595738e+00,  -1.31019838e+01,   1.52180721e+01,          3.06233294e+21],
                          [  1.81815804e+01,  -1.42908793e+01,   9.57868793e+00,         -2.84771348e+20],
                          [ -2.49671396e+01,   1.25082843e+01,  -1.71562286e+01,          2.36290426e+07],
                          [  2.67277673e+01,   1.70315414e+01,   6.12701450e+00,          7.77917232e+03],
                          [  2.49565476e+01,   2.91694684e+01,   6.29622660e+00,          2.35300027e+02],
                          [  6.11924542e+00,  -1.59943768e+00,   9.57009289e+00,          1.32906326e+11],
                          [ -1.47863653e+01,   2.41691301e+01,  -1.89981821e+01,          2.73064953e+03],
                          [  2.24070483e+01,  -2.93647433e+00,   8.19281432e+00,         -6.42000372e+17],
                          [  8.04042600e-01,   1.82710085e+01,  -1.97814534e+01,          5.48372441e-01],
                          [  1.39590390e+01,   1.97318686e+01,   2.37606635e+00,          5.51923681e+00],
                          [ -4.66640483e+00,  -2.00237930e+01,   7.40365095e+00,          4.50310752e+00],
                          [  2.76821999e+01,  -6.36563968e+00,   1.11533984e+01,         -9.28725179e+23],
                          [ -2.56764457e+01,   1.24544906e+00,   1.06407572e+01,          1.25922076e+01],
                          [  3.20447808e+00,   1.30874383e+01,   2.26098014e+01,          2.03202059e+04],
                          [ -1.24809647e+01,   4.15137113e+00,  -2.92265700e+01,          2.39621411e+08],
                          [  2.14778108e+01,  -2.35162960e+00,  -1.13758664e+01,          4.46882152e-01],
                          [ -9.85469168e+00,  -3.28157680e+00,   1.67447548e+01,         -1.07342390e+07],
                          [  1.08122310e+01,  -2.47353236e+01,  -1.15622349e+01,         -2.91733796e+03],
                          [ -2.67933347e+01,  -3.39100709e+00,   2.56006986e+01,         -5.29275382e+09],
                          [ -8.60066776e+00,  -8.02200924e+00,   1.07231926e+01,          1.33548320e+06],
                          [ -1.01724238e-01,  -1.18479709e+01,  -2.55407104e+01,          1.55436570e+00],
                          [ -3.93356771e+00,   2.11106818e+01,  -2.57598485e+01,          2.13467840e+01],
                          [  3.74750503e+00,   1.55687633e+01,  -2.92841720e+01,          1.43873509e-02],
                          [  6.99726781e+00,   2.69855571e+01,  -1.63707771e+01,          3.08098673e-02],
                          [ -2.31996011e+01,   3.47631054e+00,   9.75119815e-01,          1.79971073e-02],
                          [  2.38951044e+01,  -2.91460190e+01,  -2.50774708e+00,          9.56934814e+00],
                          [  1.52730825e+01,   5.77062507e+00,   1.21922003e+01,          1.32345307e+09],
                          [  1.74673917e+01,   1.89723426e+01,   4.94903250e+00,          9.90859484e+01],
                          [  1.88971241e+01,   2.86255413e+01,   5.52360109e-01,          1.44165360e+00],
                          [  1.02002319e+01,  -1.66855152e+01,  -2.55426235e+01,          6.56481554e+02],
                          [ -1.79474153e+01,   1.22210200e+01,  -1.84058212e+01,          8.24041812e+05],
                          [ -1.36147103e+01,   1.32365492e+00,  -7.22375200e+00,          9.92446491e+05],
                          [  7.57407832e+00,   2.59738234e+01,  -1.34139168e+01,          3.64037761e-02],
                          [  2.21110169e+00,   1.28012666e+01,   1.62529102e+01,          1.33433085e+02],
                          [ -2.64297569e+01,  -1.63176658e+01,  -1.11642006e+01,         -2.44797251e+13],
                          [ -2.46622944e+01,  -3.02147372e+00,   8.29159315e+00,         -3.21799070e+05],
                          [ -1.37215095e+01,  -1.96680183e+01,   2.91940118e+01,          3.21457520e+12],
                          [ -5.45566105e+00,   2.81292086e+01,   1.72548215e-01,          9.66973000e-01],
                          [ -1.55751298e+00,  -8.65703373e+00,   2.68622026e+01,         -3.17190834e+16],
                          [  2.45393609e+01,  -2.70571903e+01,   1.96815505e+01,          1.80708004e+37],
                          [  5.77482829e+00,   1.53203143e+01,   2.50534322e+01,          1.14304242e+06],
                          [ -1.02626819e+01,   2.36887658e+01,  -2.32152102e+01,          7.28965646e+02],
                          [ -1.30833446e+00,  -1.28310210e+01,   1.87275544e+01,         -9.33487904e+12],
                          [  5.83024676e+00,  -1.49279672e+01,   2.44957538e+01,         -7.61083070e+27],
                          [ -2.03130747e+01,   2.59641715e+01,  -2.06174328e+01,          4.54744859e+04],
                          [  1.97684551e+01,  -2.21410519e+01,  -2.26728740e+01,          3.53113026e+06],
                          [  2.73673444e+01,   2.64491725e+01,   1.57599882e+01,          1.07385118e+07],
                          [  5.73287971e+00,   1.21111904e+01,   1.33080171e+01,          2.63220467e+03],
                          [ -2.82751072e+01,   2.08605881e+01,   9.09838900e+00,         -6.60957033e-07],
                          [  1.87270691e+01,  -1.74437016e+01,   1.52413599e+01,          6.59572851e+27],
                          [  6.60681457e+00,  -2.69449855e+00,   9.78972047e+00,         -2.38587870e+12],
                          [  1.20895561e+01,  -2.51355765e+01,   2.30096101e+01,          7.58739886e+32],
                          [ -2.44682278e+01,   2.10673441e+01,  -1.36705538e+01,          4.54213550e+04],
                          [ -4.50665152e+00,   3.72292059e+00,  -4.83403707e+00,          2.68938214e+01],
                          [ -7.46540049e+00,  -1.08422222e+01,  -1.72203805e+01,         -2.09402162e+02],
                          [ -2.00307551e+01,  -7.50604431e+00,  -2.78640020e+01,          4.15985444e+19],
                          [  1.99890876e+01,   2.20677419e+01,  -2.51301778e+01,          1.23840297e-09],
                          [  2.03183823e+01,  -7.66942559e+00,   2.10340070e+01,          1.46285095e+31],
                          [ -2.90315825e+00,  -2.55785967e+01,  -9.58779316e+00,          2.65714264e-01],
                          [  2.73960829e+01,  -1.80097203e+01,  -2.03070131e+00,          2.52908999e+02],
                          [ -2.11708058e+01,  -2.70304032e+01,   2.48257944e+01,          3.09027527e+08],
                          [  2.21959758e+01,   4.00258675e+00,  -1.62853977e+01,         -9.16280090e-09],
                          [  1.61661840e+01,  -2.26845150e+01,   2.17226940e+01,         -8.24774394e+33],
                          [ -3.35030306e+00,   1.32670581e+00,   9.39711214e+00,         -1.47303163e+01],
                          [  7.23720726e+00,  -2.29763909e+01,   2.34709682e+01,         -9.20711735e+29],
                          [  2.71013568e+01,   1.61951087e+01,  -7.11388906e-01,          2.98750911e-01],
                          [  8.40057933e+00,  -7.49665220e+00,   2.95587388e+01,          6.59465635e+29],
                          [ -1.51603423e+01,   1.94032322e+01,  -7.60044357e+00,          1.05186941e+02],
                          [ -8.83788031e+00,  -2.72018313e+01,   1.88269907e+00,          1.81687019e+00],
                          [ -1.87283712e+01,   5.87479570e+00,  -1.91210203e+01,          2.52235612e+08],
                          [ -5.61338513e-01,   2.69490237e+01,   1.16660111e-01,          9.97567783e-01],
                          [ -5.44354025e+00,  -1.26721408e+01,  -4.66831036e+00,          1.06660735e-01],
                          [ -2.18846497e+00,   2.33299566e+01,   9.62564397e+00,          3.03842061e-01],
                          [  6.65661299e+00,  -2.39048713e+01,   1.04191807e+01,          4.73700451e+13],
                          [ -2.57298921e+01,  -2.60811296e+01,   2.74398110e+01,         -5.32566307e+11],
                          [ -1.11431826e+01,  -1.59420160e+01,  -1.84880553e+01,         -1.01514747e+02],
                          [  6.50301931e+00,   2.59859051e+01,  -2.33270137e+01,          1.22760500e-02],
                          [ -1.94987891e+01,  -2.62123262e+01,   3.90323225e+00,          1.71658894e+01],
                          [  7.26164601e+00,  -1.41469402e+01,   2.81499763e+01,         -2.50068329e+31],
                          [ -1.52424040e+01,   2.99719005e+01,  -2.85753678e+01,          1.31906693e+04],
                          [  5.24149291e+00,  -1.72807223e+01,   2.22129493e+01,          2.50748475e+25],
                          [  3.63207230e-01,  -9.54120862e-02,  -2.83874044e+01,          9.43854939e-01],
                          [ -2.11326457e+00,  -1.25707023e+01,   1.17172130e+00,          1.20812698e+00],
                          [  2.48513582e+00,   1.03652647e+01,  -1.84625148e+01,          6.47910997e-02],
                          [  2.65395942e+01,   2.74794672e+01,   1.29413428e+01,          2.89306132e+05],
                          [ -9.49445460e+00,   1.59930921e+01,  -1.49596331e+01,          3.27574841e+02],
                          [ -5.89173945e+00,   9.96742426e+00,   2.60318889e+01,         -3.15842908e-01],
                          [ -1.15387239e+01,  -2.21433107e+01,  -2.17686413e+01,          1.56724718e-01],
                          [ -5.30592244e+00,  -2.42752190e+01,   1.29734035e+00,          1.31985534e+00]])

        for a,b,c,expected in ref_data:
            result = hyp1f1(a,b,c)
            assert(abs(expected - result)/expected < 1e-4)

    def test_hyp1f2(self):
        pass

    def test_hyp2f0(self):
        pass

    def test_hyp2f1(self):
        # a collection of special cases taken from AMS 55
        values = [[0.5, 1, 1.5, 0.2**2, 0.5/0.2*log((1+0.2)/(1-0.2))],
                  [0.5, 1, 1.5, -0.2**2, 1./0.2*arctan(0.2)],
                  [1, 1, 2, 0.2, -1/0.2*log(1-0.2)],
                  [3, 3.5, 1.5, 0.2**2,
                      0.5/0.2/(-5)*((1+0.2)**(-5)-(1-0.2)**(-5))],
                  [-3, 3, 0.5, sin(0.2)**2, cos(2*3*0.2)],
                  [3, 4, 8, 1, gamma(8)*gamma(8-4-3)/gamma(8-3)/gamma(8-4)],
                  [3, 2, 3-2+1, -1, 1./2**3*sqrt(pi)*
                      gamma(1+3-2)/gamma(1+0.5*3-2)/gamma(0.5+0.5*3)],
                  [5, 2, 5-2+1, -1, 1./2**5*sqrt(pi)*
                      gamma(1+5-2)/gamma(1+0.5*5-2)/gamma(0.5+0.5*5)],
                  [4, 0.5+4, 1.5-2*4, -1./3, (8./9)**(-2*4)*gamma(4./3)*
                      gamma(1.5-2*4)/gamma(3./2)/gamma(4./3-2*4)],
                  # and some others
                  # ticket #424
                  [1.5, -0.5, 1.0, -10.0, 4.1300097765277476484],
                  ]
        for i, (a, b, c, x, v) in enumerate(values):
            cv = hyp2f1(a, b, c, x)
            assert_almost_equal(cv, v, 8, err_msg='test #%d' % i)

    def test_hyp3f0(self):
        pass

    def test_hyperu(self):
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

class TestBessel(TestCase):
    def test_i0(self):
        values = [[0.0, 1.0],
                  [1e-10, 1.0],
                  [0.1, 0.9071009258],
                  [0.5, 0.6450352706],
                  [1.0, 0.4657596077],
                  [2.5, 0.2700464416],
                  [5.0, 0.1835408126],
                  [20.0, 0.0897803119],
                 ]
        for i, (x, v) in enumerate(values):
            cv = i0(x) * exp(-x)
            assert_almost_equal(cv, v, 8, err_msg='test #%d' % i)

    def test_i0e(self):
        oize = i0e(.1)
        oizer = ive(0,.1)
        assert_almost_equal(oize,oizer,8)

    def test_i1(self):
        values = [[0.0, 0.0],
                  [1e-10, 0.4999999999500000e-10],
                  [0.1, 0.0452984468],
                  [0.5, 0.1564208032],
                  [1.0, 0.2079104154],
                  [5.0, 0.1639722669],
                  [20.0, 0.0875062222],
                 ]
        for i, (x, v) in enumerate(values):
            cv = i1(x) * exp(-x)
            assert_almost_equal(cv, v, 8, err_msg='test #%d' % i)

    def test_i1e(self):
        oi1e = i1e(.1)
        oi1er = ive(1,.1)
        assert_almost_equal(oi1e,oi1er,8)

    def test_iti0k0(self):
        iti0 = array(iti0k0(5))
        assert_array_almost_equal(iti0,array([31.848667776169801, 1.5673873907283657]),5)

    def test_it2i0k0(self):
        it2k = it2i0k0(.1)
        assert_array_almost_equal(it2k,array([0.0012503906973464409, 3.3309450354686687]),6)

    def test_itj0y0(self):
        it0 = array(itj0y0(.2))
        assert_array_almost_equal(it0,array([0.19933433254006822, -0.34570883800412566]),8)

    def test_it2j0y0(self):
        it2 = array(it2j0y0(.2))
        assert_array_almost_equal(it2,array([0.0049937546274601858, -0.43423067011231614]),8)

    def test_negv(self):
        assert_equal(iv(3,2), iv(-3,2))

    def test_iv(self):
        iv1 = iv(0,.1)*exp(-.1)
        assert_almost_equal(iv1,0.90710092578230106,10)

    def test_negv(self):
        assert_equal(ive(3,2), ive(-3,2))

    def test_ive(self):
        ive1 = ive(0,.1)
        iv1 = iv(0,.1)*exp(-.1)
        assert_almost_equal(ive1,iv1,10)

    def test_ivp0(self):
        assert_almost_equal(iv(1,2), ivp(0,2), 10)

    def test_ivp(self):
        y=(iv(0,2)+iv(2,2))/2
        x = ivp(1,2)
        assert_almost_equal(x,y,10)

    def test_j0(self):
        oz = j0(.1)
        ozr = jn(0,.1)
        assert_almost_equal(oz,ozr,8)

    def test_j1(self):
        o1 = j1(.1)
        o1r = jn(1,.1)
        assert_almost_equal(o1,o1r,8)

    def test_jn(self):
        jnnr = jn(1,.2)
        assert_almost_equal(jnnr,0.099500832639235995,8)

    def test_negv(self):
        assert_almost_equal(jv(-3,2), -jv(3,2), 14)

    def test_jv(self):
        values = [[0, 0.1, 0.99750156206604002],
                  [2./3, 1e-8, 0.3239028506761532e-5],
                  [2./3, 1e-10, 0.1503423854873779e-6],
                  [3.1, 1e-10, 0.1711956265409013e-32],
                  [2./3, 4.0, -0.2325440850267039],
                 ]
        for i, (v, x, y) in enumerate(values):
            yc = jv(v, x)
            assert_almost_equal(yc, y, 8, err_msg='test #%d' % i)

    def test_negv(self):
        assert_almost_equal(jve(-3,2), -jve(3,2), 14)

    def test_jve(self):
        jvexp = jve(1,.2)
        assert_almost_equal(jvexp,0.099500832639235995,8)
        jvexp1 = jve(1,.2+1j)
        z = .2+1j
        jvexpr = jv(1,z)*exp(-abs(z.imag))
        assert_almost_equal(jvexp1,jvexpr,8)

    def test_jn_zeros(self):
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

    def test_jnjnp_zeros(self):
        pass
        #jnjp = jnjnp(3)
        #assert_array_almost_equal(jnjp,(array([
        #I don't think specfun jdzo is working properly the outputs do not seem to correlate
        #to the inputs

    def test_jnp_zeros(self):
        jnp = jnp_zeros(1,5)
        assert_array_almost_equal(jnp, array([  1.84118,
                                                5.33144,
                                                8.53632,
                                                11.70600,
                                                14.86359]),4)

    def test_jnyn_zeros(self):
        jnz = jnyn_zeros(1,5)
        assert_array_almost_equal(jnz,(array([  3.83171,
                                                7.01559,
                                                10.17347,
                                                13.32369,
                                                16.47063]),
                                       array([  1.84118,
                                                5.33144,
                                                8.53632,
                                                11.70600,
                                                14.86359]),
                                       array([  2.19714,
                                                5.42968,
                                                8.59601,
                                                11.74915,
                                                14.89744]),
                                       array([  3.68302,
                                                6.94150,
                                                10.12340,
                                                13.28576,
                                                16.44006])),4)

    def test_jvp(self):
        jvprim = jvp(2,2)
        jv0 = (jv(1,2)-jv(3,2))/2
        assert_almost_equal(jvprim,jv0,10)

    def test_k0(self):
        ozk = k0(.1)
        ozkr = kv(0,.1)
        assert_almost_equal(ozk,ozkr,8)

    def test_k0e(self):
        ozke = k0e(.1)
        ozker = kve(0,.1)
        assert_almost_equal(ozke,ozker,8)

    def test_k1(self):
        o1k = k1(.1)
        o1kr = kv(1,.1)
        assert_almost_equal(o1k,o1kr,8)

    def test_k1e(self):
        o1ke = k1e(.1)
        o1ker = kve(1,.1)
        assert_almost_equal(o1ke,o1ker,8)

    def test_jacobi(self):
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

    def test_kn(self):
        kn1 = kn(0,.2)
        assert_almost_equal(kn1,1.7527038555281462,8)

    def test_negv(self):
        assert_equal(kv(3.0, 2.2), kv(-3.0, 2.2))

    def test_kv0(self):
        kv0 = kv(0,.2)
        assert_almost_equal(kv0, 1.7527038555281462, 10)
    def test_kv1(self):
        kv1 = kv(1,0.2)
        assert_almost_equal(kv1, 4.775972543220472, 10)
    def test_kv2(self):
        kv2 = kv(2,0.2)
        assert_almost_equal(kv2, 49.51242928773287, 10)


    def test_negv(self):
        assert_equal(kve(3.0, 2.2), kve(-3.0, 2.2))

    def test_kve(self):
        kve1 = kve(0,.2)
        kv1 = kv(0,.2)*exp(.2)
        assert_almost_equal(kve1,kv1,8)
        z = .2+1j
        kve2 = kve(0,z)
        kv2 = kv(0,z)*exp(z)
        assert_almost_equal(kve2,kv2,8)

    def test_kvp_v0n1(self):
        z = 2.2
        assert_almost_equal(-kv(1,z), kvp(0,z, n=1), 10)

    def test_kvp_n1(self):
        v = 3.
        z = 2.2
        xc = -kv(v+1,z) + v/z*kv(v,z)
        x = kvp(v,z, n=1)
        assert_almost_equal(xc, x, 10)   #this function (kvp) is broken

    def test_kvp_n2(self):
        v = 3.
        z = 2.2
        xc = (z**2+v**2-v)/z**2 * kv(v,z) + kv(v+1,z)/z
        x = kvp(v, z, n=2)
        assert_almost_equal(xc, x, 10)

    def test_y0(self):
        oz = y0(.1)
        ozr = yn(0,.1)
        assert_almost_equal(oz,ozr,8)

    def test_y1(self):
        o1 = y1(.1)
        o1r = yn(1,.1)
        assert_almost_equal(o1,o1r,8)

    def test_y0_zeros(self):
        yo,ypo = y0_zeros(2)
        zo,zpo = y0_zeros(2,complex=1)
        all = r_[yo,zo]
        allval = r_[ypo,zpo]
        assert_array_almost_equal(abs(yv(0.0,all)),0.0,11)
        assert_array_almost_equal(abs(yv(1,all)-allval),0.0,11)


    def test_y1_zeros(self):
        y1 = y1_zeros(1)
        assert_array_almost_equal(y1,(array([2.19714]),array([0.52079])),5)

    def test_y1p_zeros(self):
        y1p = y1p_zeros(1,complex=1)
        assert_array_almost_equal(y1p,(array([ 0.5768+0.904j]), array([-0.7635+0.5892j])),3)

    def test_yn_zeros(self):
        an = yn_zeros(4,2)
        assert_array_almost_equal(an,array([ 5.64515,  9.36162]),5)

    def test_ynp_zeros(self):
        ao = ynp_zeros(0,2)
        assert_array_almost_equal(ao,array([ 2.19714133, 5.42968104]),6)

    def test_yn(self):
        yn2n = yn(1,.2)
        assert_almost_equal(yn2n,-3.3238249881118471,8)

    def test_negv(self):
        assert_almost_equal(yv(-3,2), -yv(3,2), 14)

    def test_yv(self):
        yv2 = yv(1,.2)
        assert_almost_equal(yv2,-3.3238249881118471,8)

    def test_negv(self):
        assert_almost_equal(yve(-3,2), -yve(3,2), 14)

    def test_yve(self):
        yve2 = yve(1,.2)
        assert_almost_equal(yve2,-3.3238249881118471,8)
        yve2r = yv(1,.2+1j)*exp(-1)
        yve22 = yve(1,.2+1j)
        assert_almost_equal(yve22,yve2r,8)

    def test_yvp(self):
        yvpr = (yv(1,.2) - yv(3,.2))/2.0
        yvp1 = yvp(2,.2)
        assert_array_almost_equal(yvp1,yvpr,10)


class TestLaguerre(TestCase):
    def test_laguerre(self):
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

    def test_genlaguerre(self):
        k = 5*rand()-0.9
        lag0 = genlaguerre(0,k)
        lag1 = genlaguerre(1,k)
        lag2 = genlaguerre(2,k)
        lag3 = genlaguerre(3,k)
        assert_equal(lag0.c,[1])
        assert_equal(lag1.c,[-1,k+1])
        assert_almost_equal(lag2.c,array([1,-2*(k+2),(k+1.)*(k+2.)])/2.0)
        assert_almost_equal(lag3.c,array([-1,3*(k+3),-3*(k+2)*(k+3),(k+1)*(k+2)*(k+3)])/6.0)


# Base polynomials come from Abrahmowitz and Stegan
class TestLegendre(TestCase):
    def test_legendre(self):
        leg0 = legendre(0)
        leg1 = legendre(1)
        leg2 = legendre(2)
        leg3 = legendre(3)
        leg4 = legendre(4)
        leg5 = legendre(5)
        assert_equal(leg0.c,[1])
        assert_equal(leg1.c,[1,0])
        assert_equal(leg2.c,array([3,0,-1])/2.0)
        assert_almost_equal(leg3.c,array([5,0,-3,0])/2.0)
        assert_almost_equal(leg4.c,array([35,0,-30,0,3])/8.0)
        assert_almost_equal(leg5.c,array([63,0,-70,0,15,0])/8.0)


class TestLambda(TestCase):
    def test_lmbda(self):
        lam = lmbda(1,.1)
        lamr = (array([jn(0,.1), 2*jn(1,.1)/.1]),
                array([jvp(0,.1), -2*jv(1,.1)/.01 + 2*jvp(1,.1)/.1]))
        assert_array_almost_equal(lam,lamr,8)

class TestLog1p(TestCase):
    def test_log1p(self):
        l1p = (log1p(10),log1p(11),log1p(12))
        l1prl = (log(11),log(12),log(13))
        assert_array_almost_equal(l1p,l1prl,8)

    def test_log1pmore(self):
        l1pm = (log1p(1),log1p(1.1),log1p(1.2))
        l1pmrl = (log(2),log(2.1),log(2.2))
        assert_array_almost_equal(l1pm,l1pmrl,8)

class TestLegendreFunctions(TestCase):
    def test_lpmn(self):
        lp = lpmn(0,2,.5)
        assert_array_almost_equal(lp,(array([       [ 1.00000 ,
                                                      0.50000,
                                                      -0.12500]]),
                                      array([       [ 0.00000 ,
                                                      1.00000 ,
                                                      1.50000]])),4)

    def test_lpn(self):
        lpnf = lpn(2,.5)
        assert_array_almost_equal(lpnf,(array(      [ 1.00000 ,
                                                        0.50000,
                                                        -0.12500]),
                                      array(        [ 0.00000 ,
                                                      1.00000 ,
                                                      1.50000])),4)

    def test_lpmv(self):
        lp = lpmv(0,2,.5)
        assert_almost_equal(lp,-0.125,3)

    def test_lqmn(self):
        lqmnf = lqmn(0,2,.5)
        lqmnf = lqmn(0,2,.5)
        lqf = lqn(2,.5)
        assert_array_almost_equal(lqmnf[0][0],lqf[0],4)
        assert_array_almost_equal(lqmnf[1][0],lqf[1],4)

    def test_lqn(self):
        lqf = lqn(2,.5)
        assert_array_almost_equal(lqf,(array([ 0.5493, -0.7253, -0.8187]),
                                       array([ 1.3333,  1.216 , -0.8427])),4)

class TestMathieu(TestCase):

    def test_mathieu_a(self):
        pass

    def test_mathieu_even_coef(self):
        mc =  mathieu_even_coef(2,5)
        #Q not defined broken and cannot figure out proper reporting order

    def test_mathieu_odd_coef(self):
        pass
            #same problem as above

class TestFresnelIntegral(TestCase):

    def test_modfresnelp(self):
        pass

    def test_modfresnelm(self):
        pass

class TestOblCvSeq(TestCase):
    def test_obl_cv_seq(self):
        obl = obl_cv_seq(0,3,1)
        assert_array_almost_equal(obl,array([ -0.348602,
                                              1.393206,
                                              5.486800,
                                              11.492120]),5)

class TestParabolicCylinder(TestCase):
    def test_pbdn_seq(self):
        pb = pbdn_seq(1,.1)
        assert_array_almost_equal(pb,(array([ 0.9975,
                                              0.0998]),
                                      array([-0.0499,
                                             0.9925])),4)

    def test_pbdv(self):
        pbv = pbdv(1,.2)
        derrl = 1/2*(.2)*pbdv(1,.2)[0] - pbdv(0,.2)[0]

    def test_pbdv_seq(self):
        pbn = pbdn_seq(1,.1)
        pbv = pbdv_seq(1,.1)
        assert_array_almost_equal(pbv,(real(pbn[0]),real(pbn[1])),4)

class TestPolygamma(TestCase):
    # from Table 6.2 (pg. 271) of A&S
    def test_polygamma(self):
        poly2 = polygamma(2,1)
        poly3 = polygamma(3,1)
        assert_almost_equal(poly2,-2.4041138063,10)
        assert_almost_equal(poly3,6.4939394023,10)

class TestProCvSeq(TestCase):
    def test_pro_cv_seq(self):
        prol = pro_cv_seq(0,3,1)
        assert_array_almost_equal(prol,array([  0.319000,
                                               2.593084,
                                               6.533471,
                                               12.514462]),5)

class TestPsi(TestCase):
    def test_psi(self):
        ps = psi(1)
        assert_almost_equal(ps,-0.57721566490153287,8)

class TestRadian(TestCase):
    def test_radian(self):
        rad = radian(90,0,0)
        assert_almost_equal(rad,pi/2.0,5)

    def test_radianmore(self):
        rad1 = radian(90,1,60)
        assert_almost_equal(rad1,pi/2+0.0005816135199345904,5)

class TestRiccati(TestCase):
    def test_riccati_jn(self):
        jnrl = (sph_jn(1,.2)[0]*.2,sph_jn(1,.2)[0]+sph_jn(1,.2)[1]*.2)
        ricjn = riccati_jn(1,.2)
        assert_array_almost_equal(ricjn,jnrl,8)

    def test_riccati_yn(self):
        ynrl = (sph_yn(1,.2)[0]*.2,sph_yn(1,.2)[0]+sph_yn(1,.2)[1]*.2)
        ricyn = riccati_yn(1,.2)
        assert_array_almost_equal(ricyn,ynrl,8)

class TestRound(TestCase):
    def test_round(self):
        rnd = map(int,(round(10.1),round(10.4),round(10.5),round(10.6)))

        # Note: According to the documentation, scipy.special.round is
        # supposed to round to the nearest even number if the fractional
        # part is exactly 0.5. On some platforms, this does not appear
        # to work and thus this test may fail. However, this unit test is
        # correctly written.
        rndrl = (10,10,10,11)
        assert_array_equal(rnd,rndrl)

class _test_sh_legendre(TestCase):

    def test_sh_legendre(self):
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

class _test_sh_chebyt(TestCase):

    def test_sh_chebyt(self):
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


class _test_sh_chebyu(TestCase):

    def test_sh_chebyu(self):
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

class _test_sh_jacobi(TestCase):

    def test_sh_jacobi(self):
        # G^(p,q)_n(x) = n! gamma(n+p)/gamma(2*n+p) * P^(p-q,q-1)_n(2*x-1)
        conv = lambda n,p: _gam(n+1)*_gam(n+p)/_gam(2*n+p)
        psub = poly1d([2,-1])
        q = 4*rand()
        p = q-1 + 2*rand()
        #print "shifted jacobi p,q = ", p, q
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

class TestSpherical(TestCase):
    def test_sph_harm(self):
        pass

    def test_sph_in(self):
        i1n = sph_in(1,.2)
        inp0 = (i1n[0][1])
        inp1 = (i1n[0][0] - 2.0/0.2 * i1n[0][1])
        assert_array_almost_equal(i1n[0],array([1.0066800127054699381,
                                                0.066933714568029540839]),12)
        assert_array_almost_equal(i1n[1],[inp0,inp1],12)

    def test_sph_inkn(self):
        spikn = r_[sph_in(1,.2)+sph_kn(1,.2)]
        inkn = r_[sph_inkn(1,.2)]
        assert_array_almost_equal(inkn,spikn,10)

    def test_sph_jn(self):
        s1 = sph_jn(2,.2)
        s10 = -s1[0][1]
        s11 = s1[0][0]-2.0/0.2*s1[0][1]
        s12 = s1[0][1]-3.0/0.2*s1[0][2]
        assert_array_almost_equal(s1[0],[0.99334665397530607731,
                                      0.066400380670322230863,
                                      0.0026590560795273856680],12)
        assert_array_almost_equal(s1[1],[s10,s11,s12],12)

    def test_sph_jnyn(self):
        jnyn = r_[sph_jn(1,.2) + sph_yn(1,.2)]  # tuple addition
        jnyn1 = r_[sph_jnyn(1,.2)]
        assert_array_almost_equal(jnyn1,jnyn,9)

    def test_sph_kn(self):
        kn = sph_kn(2,.2)
        kn0 = -kn[0][1]
        kn1 = -kn[0][0]-2.0/0.2*kn[0][1]
        kn2 = -kn[0][1]-3.0/0.2*kn[0][2]
        assert_array_almost_equal(kn[0],[6.4302962978445670140,
                                         38.581777787067402086,
                                         585.15696310385559829],12)
        assert_array_almost_equal(kn[1],[kn0,kn1,kn2],9)

    def test_sph_yn(self):
        sy1 = sph_yn(2,.2)[0][2]
        sy2 = sph_yn(0,.2)[0][0]
        sphpy = (sph_yn(1,.2)[0][0]-2*sph_yn(2,.2)[0][2])/3 #correct derivative value
        assert_almost_equal(sy1,-377.52483,5)#previous values in the system
        assert_almost_equal(sy2,-4.9003329,5)
        sy3 = sph_yn(1,.2)[1][1]
        assert_almost_equal(sy3,sphpy,4) #compare correct derivative val. (correct =-system val).

if __name__ == "__main__":
    nose.run(argv=['', __file__])
