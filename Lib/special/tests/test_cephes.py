#! /usr/bin/env python 
#

import sys
from scipy_test.testing import set_package_path
set_package_path()
from special import cephes
del sys.path[0]

import unittest
from scipy_test.testing import assert_array_equal, assert_array_almost_equal
from scipy_test.testing import assert_almost_equal,assert_equal
from scipy_test.testing import ScipyTestCase

class test_cephes(ScipyTestCase):
    def check_airy(self):
        print "***********Yes, this is testing******************"
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
        assert_equal(cephes.cbrt(1),1.0)

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
        assert_equal(cephes.exp10(2),100.0)
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
        assert_equal(cephes.nrdtrimn(0.5,1,1),1.0)
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
        assert_equal(cephes.pro_ang1_cv(1,1,1,1,0),(1.0,0.0))
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

#-------------------------------------------
def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append(unittest.makeSuite(test_cephes,'check_'))
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner

if __name__ == "__main__":
    if len(sys.argv)>1:
        level = eval(sys.argv[1])
    else:
        level = 1
    test(level)


#Old contents
'''
import sys
sys.path.append('..')

import cephes,Test

import glob, re

if sys.argv[-1] == '-b': reference=(0==0)
else: reference=(0==1)

for filenam in glob.glob('fncs_*.dat'):
    file=open(filenam,'r')
    _last=0
    while 1:
        testline=""
        while 1:
            line=file.readline()
            if line=="":
                _last=1
                break
            if line[-1:] == '\n': line = line[:-1]
            if line[-1:] == '\r': line = line[:-1]
            line=re.sub("#.*","",line)
            if re.search('\S',line):
                testline=testline+line
            else: break
        if testline != "":
            a=eval('Test.Test('+testline+',ref=reference)')
            rel_err=a.test()
            print a.name,
            if rel_err>0: print 'MAX_ERR: ', rel_err*100,'%'
            else: print 'PASS'
        if _last: break
    
'''
