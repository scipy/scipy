#this program corresponds to special.py
import os
import sys
#import fpformat
import unittest
import scipy_base.limits as limits
from scipy_test.testing import assert_array_equal, assert_equal
from scipy_test.testing import assert_almost_equal, assert_array_almost_equal
from scipy_base import *

from scipy_test.testing import set_package_path
set_package_path()
from special import *
del sys.path[0]


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
					  0.51053920]),8)
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
        for n in range(4):
            b1[n] = b[n]*exp(2.0/3.0*0.01*sqrt(0.01))
        assert_array_almost_equal(a,b1,5)

class test_arange(unittest.TestCase):

    def check_arange(self):
        numstring = arange(0,2.2,.1)
        assert_equal(numstring,array([0.,0.1,0.2,0.3,
                                      0.4,0.5,0.6,0.7,
                                      0.8,0.9,1.,1.1,
                                      1.2,1.3,1.4,1.5,
                                      1.6,1.7,1.8,1.9,
                                      2.,2.1]))
        numstringa = arange(3,4,.3)
        assert_array_equal(numstringa, array([3.,3.3,3.6,3.9]))
        numstringb = arange(3,27,3)
        assert_array_equal(numstringb,array([3.,6.,9.,12.,
                                             15.,18.,21.,24.,27.]))
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
        assert_array_equal(a,z)

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
        assert_array_almost_equal(bi,array([5.02622,
                                            9.45541,
                                            13.89349,
                                            18.33398,
                                            22.77544]),4)

class test_beip_zeros(unittest.TestCase):

    def check_beip_zeros(self):
        bip = beip_zeros(5)
        assert_array_almost_equal(bip,array([3.77320,
                                             8.28099,
                                             12.74215,
                                             17.19343,
                                             21.64114]),4)

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
        chebc = chebyc(1)(.2)
        assert_almost_equal(chebc,.2000,4)

class test_chebys(unittest.TestCase):

    def check_chebys(self):
        chebs = chebys(1)(.2)
        assert_almost_equal(chebs,.2000,4)

class test_chebyt(unittest.TestCase):

    def check_chebyt(self):
        cheb = chebyt(1)(.2)
        assert_equal(cheb,0.2)

class test_chebyu(unittest.TestCase):

    def check_chebyu(self):
        chebu = chebyu(1)(.2)
        assert_equal(chebu,0.4)

class test_choose(unittest.TestCase):

    def check_choose(self):
        carray = [1,3,2,4,6,5]
        chose = choose([1,3,5],carray)
        assert_array_equal(chose,array([3,4,5]))

class test_cbrt(unittest.TestCase):

    def check_cbrt(self):
        cb = cbrt(27)
        cbrl = 27**(1.0/3.0)
        assert_equal(cb,cbrt)

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
        el = array(ellipj(0.2))
        rel = array([sin(0.2),cos(0.2),1.0,0.20])
        assert_array_almost_equal(el,rel,6)

class test_ellipk(unittest.TestCase):

    def check_ellipk(self):
        elk = ellipk(.2)
        assert_almost_equal(elk,0.1659623598610528,8)

class test_ellipkinc(unittest.TestCase):

    def check_ellipkinc(self):
        elkinc = ellipkinc(pi/2,.2)
        assert_almost_equal(elkinc,2.2572053268208538,8)

class test_ellipe(unittest.TestCase):

    def check_ellipe(self):
        ele = ellipe(.2)
        assert_almost_equal(ele,1.4890350580958529,8)

class test_ellipinc(unittest.TestCase):

    def check_ellipinc(self):
        eleinc = ellipeinc(pi/2,.2)
        assert_almost_equal(eleinc,1.1784899243278384,8)

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
        eu1 = euler(5)
        assert_array_equal(eu1,array([1,0,-1,0,5,0,-61]))
        eu2 = euler(20)
        assert_almost_equal(eu2[20]/1e14,3.70371188237548,12)

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
        assert_equal(ex,exrl)

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

    def check_fresnel_zeros(self):
        fr = fresnel_zeros(4)
        assert_array_almost_equal(fr,(array([ 2.0093+0.2886j,
                                             2.8335+0.2443j,
                                             3.4675+0.2185j,
                                             4.0026+0.2008j]),
                                  array([ 1.7437+0.3057j,
                                          2.6515+0.2529j,
                                          3.3208+0.2239j ,
                                          3.8759+0.2047j])),4)

class test_fresnelc_zeros(unittest.TestCase):

    def check_fresnelc_zeros(self):
        frc = fresnelc_zeros(4)
        asser_array_almost_equal(frc,array([ 1.7437+0.3057j,
                                          2.6515+0.2529j,
                                          3.3208+0.2239j ,
                                          3.8759+0.2047j]),4)

class test_fresnels_zeros(unittest.TestCase):

    def check_fresnels_zeros(self):
        frs = fresnels_zeros(4)
        assert_array_almost_equal(frs,array([ 2.0093+0.2886j,
                                             2.8335+0.2443j,
                                             3.4675+0.2185j,
                                             4.0026+0.2008j]),4)


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

class test_genlaguerre(unittest.TestCase):

    def check_genlaguerre(self):
        gen = genlaguerre(1,0)(.5)
        assert_equal(gen,0.5)

class test_hankel1(unittest.TestCase):

    def check_hankel1(self):
        hank1 = hankel1(1,.1)
        hankrl = (jv(1,.1)+yv(1,.1)*1j)
        assert_almost_equal(hank1,hankrl,8)

class test_hankel1e(unittest.TestCase):

    def check_hankel1e(self):
        hank1e = hankel1e(1,.1)
        hankrle = hankel1(1,.1)*exp(.1j)
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

        herm = hermite(3)(.5)
        assert_almost_equal(herm,-5.0,8)

class test_gegenbauer(unittest.TestCase):

    def check_gegenbauer(self):
        geg = gegenbauer(2,.5)(0)
        assert_equal(geg,-0.5)

class test_h1vp(unittest.TestCase):

    def check_h1vp(self):
        
        h1 = h1vp(1,.1)
        h1real = (jvp(1,.1)+yvp(1,.1)*1j)
        assert_almost_equal(h1,h1real,8)

class test_h2vp(unittest.TestCase):

    def check_h2vp(self):
        h2 = h2vp(1,.1)
        h2real = (jvp(1,.1)+yvp(1,.1)*1j)
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
        hypu = hyperu(2,2,.5)
        hprl = (pi/sin(pi/b))*(hyp1f1(2,2,.5)/(gamma(1)*gamma(2))-.5**(-1)*hyp1f1(1,0,.5)/(gamma(2)))
        assert_almost_equal(hypu,hprl,8)

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
        y=(iv(0,2)+iv(2,2))/2
        x = ivp(1,2,1)
        assert_almost_equal(x,y,4)

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
        j = jacobi(1,1,1)(.2)

class test_jn(unittest.TestCase):

    def check_jn(self):
        jnnr = jn(1,.2)
        assert_almost_equal(jnnr,0.099500832639235995,8)

class test_jv(unittest.TestCase):

    def check_jv(self):
        jc = jv(0,.1)
        assert_almost_equal(jnnr,0.99750156206604002,8)

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
        jv = (jv(1,2)-jv(3,2))/2
        assert_almost_equal(jvprim,jv,4)

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

    def check_kelvin_zeros(self):
        kelvz = kelvin_zeros(5)
        assert_array_almost_equal(kelvz,(array([ 2.84892,
                                                 7.23883,
                                                 11.67396,
                                                 16.11356,
                                                 20.55463]),
                                        array([	 5.02622,
                                                 9.45541,
                                                 13.89349,
                                                 18.33398,
                                                 22.77544]),
                                        array([	 1.71854,
                                                 6.12728,
                                                 10.56294,
                                                 15.00269,
                                                 19.44381]),
                                        array([	 3.91467,
                                                 8.34422,
                                                 12.78256,
                                                 17.22314,
                                                 21.66464]),
                                        array([	 6.03871,
                                                 10.51364,
                                                 14.96844,
                                                 19.41758,
                                                 23.86430]),
                                        array([	 3.77320,
                                                 8.28099,
                                                 12.74215,
                                                 17.19343,
                                                 21.64114]),
                                        array([	 2.66584,
                                                 7.17212,
                                                 11.63218,
                                                 16.08312,
                                                 20.53068]),
                                        array([	 4.93181,
                                                 9.40405,
                                                 13.85827,
                                                 18.30717,
                                                 22.75379])),4)

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
        assert_array_almost_equal(ker,array([  2.66584,
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
        kvprimrl = (kv(0,2) + kv(2,2))/2
        assert_almost_equal(kvprim,kvprimrl,4)	 #this function (kvp) is broken

class test_laguerre(unittest.TestCase):

    def check_laguerre(self):
        lag = laguerre(1)(.5)
        assert_equal(lag,.5)

class test_legendre(unittest.TestCase):

    def check_legendre(self):
        a = legendre(2)(2)
        assert_equal(a,5.5)

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
        assert_array_almost_equal(lpnf,(array([	      [ 1.00000 ,
                                                        0.50000,
                                                        -0.12500]]),
                                      array([	    [ 0.00000 ,
                                                      1.00000 ,
                                                      1.50000]])),4)

class test_lpmv(unittest.TestCase):

    def check_lpmv(self):
        lp = lpmv(0,2,.5)
        assert_almost_equal(lp,-0.125,3)

class test_lqmn(unittest.TestCase):

    def check_lqmn(self):
        disp("here")        
        lqmnf = lqmn(0,2,.5)
        lqmnf = lqmn(0,2,.5)
        lqf = lqn(2,.5)
        disp("done")
        assert_array_almost_equal(lqmnf,lqf,4)

                


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
                                              11.492120]),6)

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

class test_pbdv_seq(unittest.TestCase):

    def check_pbdv_seq(self):
        pbn = pbdn_seq(1,.1)
        pbv = pbdv_seq(1,.1)
        assert_array_almost_equal(pbv,(real(pbn[0]),real(pbn[1])),4)

class test_pbvv_seq(unittest.TestCase):

    def check_pbvv_seq(self):
        pass

class test_polygamma(unittest.TestCase):

    def check_polygamma(self):
        scipy_base.disp("Here.")
        poly2 = polygamma(2,1)
        poly3 = polygamma(3,1)
        assert_almost_equal(poly2,-2.4041138063191885,10)
        assert_almost_equal(poly2,6.4939394022668289,10)

class test_pro_cv_seq(unittest.TestCase):

    def check_pro_cv_seq(self):
        prol = pro_cv_seq(0,3,1)
        assert_array_almost_equal(prol,array([  0.319000,
                                               2.593084,
                                               6.533471,
                                               12.514462]),6)

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
        b = reshape(a,(3,2))
        assert_array_equal(b,array([[1, 2, 3],
                                    [4, 5, 6]]))
        c = reshape(a,(3,2))
        assert_array_equal(b,array([[1, 2],
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
        ricyn = riccati_jn(1,.2)
        assert_array_almost_equal(ricyn,ynrl,8)

class test_round(unittest.TestCase):

    def check_round(self):
        rnd = (round(10.1),round(10.4),round(10.5),round(10.6))
        rndrl = (10,10,11,11)
        assert_array_equal(rnd,rndrl)

class test_sh_legendre(unittest.TestCase):

    def check_sh_legendre(self):
        shleg = sh_legendre(1)(4)
        re = legendre(1)(4)-0.5
        assert_almost_equal(shleg,re,4)

class test_sh_chebyt(unittest.TestCase):

    def check_sh_chebyt(self):
        shbyt = sh_chebyt(1)(1)
        shre = chebyt(1)(1)-0.5
        assert_almost_equal(shbyt,shre,4)

class test_sh_chebyu(unittest.TestCase):

    def check_sh_chebyu(self):
        shbyu = sh_chebyu(1)(.5)
        chbre = chebyu(1)(.5)-0.5
        assert_almost_equal(shbyu,chbre,4)

class test_sh_jacobi(unittest.TestCase):

    def check_sh_jacobi(self):
        shjc = sh_jacobi(1,1,1)(.2)
        jcre = jacobi(1,1,1)(.2)-0.5
        assert_array_almost_equal(shjc,jcre,4)

class test_sinc(unittest.TestCase):

    def check_sinc(self):
        c = arange(-2,2,.1)
        y = sinc(c)
        yre = sin(pi*c)/(pi*c)
	yre[20] = 1.0
        assert_array_almost_equal(y, yre,4)

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
        i1n = sph_kn(1,.2)
        inp = (sph_kn(1,.2)[0][0]+2*sph_kn(2,.2)[0][2])/3
        assert_array_almost_equal(i1n[0],array([1.00668001,0.06693370]),8)
        assert_almost_equal(i1n[1][1],inp)

class test_sph_inkn(unittest.TestCase):

    def check_sph_inkn(self):
        spikn = (sph_in(1,.2),sph_kn(1,.2))
        inkn = sph_inkn(1,.2)
        assert_array_almost_equal(inkn,spikn,8)

class test_sph_jn(unittest.TestCase):

    def check_sph_jn(self):
        s1 = sph_jn(2,.2)[0][2]
        s2 = sph_jn(0,.2)[0][0]
        sphp = (sph_jn(1,.2)[0][0]-2*sph_jn(2,.2)[0][2])/3
        assert_almost_equal(s1,0.0059615249,8)
        assert_almost_equal(s2,0.99334665,8)
        s3 = sph_jn(1,.2)[1][1]
        assert_almost_equal(s3,sphp,4)

class test_sph_jnyn(unittest.TestCase):

    def check_sph_jnyn(self):
        jnyn = (sph_jn(1,.2),sph_jn(1,.2))
        jnyn1 = sph_jnyn(1,.2)
        assert_almost_equal(jnyn1,jnyn,8)

class test_sph_kn(unittest.TestCase):

    def check_sph_kn(self):
        kn = sph_kn(1,.2)
        knp = (sph_kn(1,.2)[0][0]+2*sph_kn(2,.2)[0][2])/3
        assert_array_almost_equal(kn[0],array([	 6.43029630,  38.5817778]),8)
        assert_almost_equal(kn[1][1],knp)

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
        yo = y0_zeros(1)
        assert_array_almost_equal(yo,array(array([0.89357697+0.j]),array([ 0.87942080+0.j])),5)

class test_y1_zeros(unittest.TestCase):

    def check_y1_zeros(self):
        y1 = y1_zeros(1)
        assert_array_almost_equal(y1,(array([2.19714]),array([0.52079])),5)
        
class test_y1p_zeros(unittest.TestCase):

    def check_y1p_zeros(self):
        y1p = y1p_zeros(1)
        assert_array_almost_equal(y1p,(array([ 0.5768+0.904j]), array([-0.7635+0.5892j])),5)

class test_yn_zeros(unittest.TestCase):

    def check_yn_zeros(self):
        an = yn_zeros(4,2)
        assert_array_almost_equal(an,array([ 5.64515,  9.36162]),5)

class test_ynp_zeros(unittest.TestCase):

    def check_ynp_zeros(self):
        ao = ynp_zeros(0,2)
        assert_array_almost_equal(ao,array([ 0.87942080, -0.40254267]),6)

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

class test_general_function(unittest.TestCase):

    def check_simple(self):
        def addsubtract(a,b):
            if a > b:
                return a - b
            else:
                return a + b
        f = general_function(addsubtract)
        r = f([0,3,6,9],[1,3,5,7])
        assert_array_equal(r,[1,6,1,2])

####### Testing ##############


### Means test is not done yet
#E   Means test is giving error (E)
#F   Means test is failing (F)
#EF  Means test is giving error and Failing
#!   Means test is segfaulting

def test_suite(level=1):
    suites = []
    if level > 0:
#F	suites.append( unittest.makeSuite(test_airy,'check_') )
#F	suites.append( unittest.makeSuite(test_airye,'check_') )
#F	suites.append( unittest.makeSuite(test_arange,'check_') )
	suites.append( unittest.makeSuite(test_ai_zeros,'check_') )
#F	suites.append( unittest.makeSuite(test_array,'check_') )
	suites.append( unittest.makeSuite(test_assoc_laguerre,'check_') )
###	suites.append( unittest.makeSuite(test_besselpoly,'check_') )
	suites.append( unittest.makeSuite(test_bei,'check_') )
	suites.append( unittest.makeSuite(test_beip,'check_') )
	suites.append( unittest.makeSuite(test_ber,'check_') )
	suites.append( unittest.makeSuite(test_berp,'check_') )
#F	suites.append( unittest.makeSuite(test_bei_zeros,'check_') )
#F	suites.append( unittest.makeSuite(test_beip_zeros,'check_') )
	suites.append( unittest.makeSuite(test_ber_zeros,'check_') )
	suites.append( unittest.makeSuite(test_berp_zeros,'check_') )
	suites.append( unittest.makeSuite(test_bernoulli,'check_') )
	suites.append( unittest.makeSuite(test_beta,'check_') )
	suites.append( unittest.makeSuite(test_betaln,'check_') )
	suites.append( unittest.makeSuite(test_betainc,'check_') )
	suites.append( unittest.makeSuite(test_betaincinv,'check_') )
	suites.append( unittest.makeSuite(test_bi_zeros,'check_') )
	suites.append( unittest.makeSuite(test_chebyc,'check_') )
	suites.append( unittest.makeSuite(test_chebys,'check_') )
#F	suites.append( unittest.makeSuite(test_chebyt,'check_') )
#F	suites.append( unittest.makeSuite(test_chebyu,'check_') )
	suites.append( unittest.makeSuite(test_choose,'check_') )
#F	suites.append( unittest.makeSuite(test_cbrt,'check_') )
	suites.append( unittest.makeSuite(test_cosdg,'check_') )
	suites.append( unittest.makeSuite(test_cosm1,'check_') )
	suites.append( unittest.makeSuite(test_cotdg,'check_') )
#E	suites.append( unittest.makeSuite(test_ellipj,'check_') )	 
#F      suites.append( unittest.makeSuite(test_ellipk,'check_') )
#F	suites.append( unittest.makeSuite(test_ellipkinc,'check_') )
	suites.append( unittest.makeSuite(test_ellipe,'check_') )
#F	suites.append( unittest.makeSuite(test_ellipinc,'check_') )
	suites.append( unittest.makeSuite(test_erf,'check_') )
	suites.append( unittest.makeSuite(test_erf_zeros,'check_') )
	suites.append( unittest.makeSuite(test_erfcinv,'check_') )
	suites.append( unittest.makeSuite(test_errprint,'check_') )
#F	suites.append( unittest.makeSuite(test_euler,'check_') )
	suites.append( unittest.makeSuite(test_exp2,'check_') )	       
	suites.append( unittest.makeSuite(test_exp10,'check_') )
	suites.append( unittest.makeSuite(test_expm1,'check_') )
	suites.append( unittest.makeSuite(test_fresnel,'check_') )
#F	suites.append( unittest.makeSuite(test_fresnel_zeros,'check_') )
#E	suites.append( unittest.makeSuite(test_fresnelc_zeros,'check_') )
#F	suites.append( unittest.makeSuite(test_fresnels_zeros,'check_') )
	suites.append( unittest.makeSuite(test_gamma,'check_') )
	suites.append( unittest.makeSuite(test_gammaln,'check_') )
	suites.append( unittest.makeSuite(test_gammainc,'check_') )
        suites.append( unittest.makeSuite(test_gammaincc,'check_') )
	suites.append( unittest.makeSuite(test_gammaincinv,'check_') )
        suites.append( unittest.makeSuite(test_gammainccinv,'check_') )	      
	suites.append( unittest.makeSuite(test_genlaguerre,'check_') )
	suites.append( unittest.makeSuite(test_hankel1,'check_') )
#F	suites.append( unittest.makeSuite(test_hankel1e,'check_') )
	suites.append( unittest.makeSuite(test_hankel2,'check_') )
	suites.append( unittest.makeSuite(test_hankel2e,'check_') )
	suites.append( unittest.makeSuite(test_hermite,'check_') )
	suites.append( unittest.makeSuite(test_gegenbauer,'check_') )
	suites.append( unittest.makeSuite(test_h1vp,'check_') )
#F	suites.append( unittest.makeSuite(test_h2vp,'check_') )
	suites.append( unittest.makeSuite(test_hyp0f1,'check_') )	 
	suites.append( unittest.makeSuite(test_hyp1f2,'check_') )
	suites.append( unittest.makeSuite(test_hyp2f0,'check_') )
	suites.append( unittest.makeSuite(test_hyp2f1,'check_') )
	suites.append( unittest.makeSuite(test_hyp3f0,'check_') )
#E	suites.append( unittest.makeSuite(test_hyperu,'check_') )
	suites.append( unittest.makeSuite(test_i0,'check_') )
	suites.append( unittest.makeSuite(test_i0e,'check_') )
	suites.append( unittest.makeSuite(test_i1,'check_') )
	suites.append( unittest.makeSuite(test_i1e,'check_') )
	suites.append( unittest.makeSuite(test_iti0k0,'check_') )
	suites.append( unittest.makeSuite(test_it2i0k0,'check_') )	  
	suites.append( unittest.makeSuite(test_itj0y0,'check_') )
	suites.append( unittest.makeSuite(test_it2j0y0,'check_') )
	suites.append( unittest.makeSuite(test_iv,'check_') )
	suites.append( unittest.makeSuite(test_ive,'check_') )
#F	suites.append( unittest.makeSuite(test_ivp,'check_') )
	suites.append( unittest.makeSuite(test_j0,'check_') )
	suites.append( unittest.makeSuite(test_j1,'check_') )
	suites.append( unittest.makeSuite(test_jacobi,'check_') )
	suites.append( unittest.makeSuite(test_jn,'check_') )
#E	suites.append( unittest.makeSuite(test_jv,'check_') )	     
	suites.append( unittest.makeSuite(test_jve,'check_') )
	suites.append( unittest.makeSuite(test_jn_zeros,'check_') )
###	 suites.append( unittest.makeSuite(test_jnjnp_zeros,'check_') )
	suites.append( unittest.makeSuite(test_jnp_zeros,'check_') )
	suites.append( unittest.makeSuite(test_jnyn_zeros,'check_') )
#E	suites.append( unittest.makeSuite(test_jvp,'check_') )
	suites.append( unittest.makeSuite(test_k0,'check_') )
	suites.append( unittest.makeSuite(test_k0e,'check_') )
	suites.append( unittest.makeSuite(test_k1,'check_') )
	suites.append( unittest.makeSuite(test_k1e,'check_') )	      
	suites.append( unittest.makeSuite(test_kei,'check_') )
	suites.append( unittest.makeSuite(test_kelvin,'check_') )
	suites.append( unittest.makeSuite(test_keip,'check_') )
	suites.append( unittest.makeSuite(test_ker,'check_') )
	suites.append( unittest.makeSuite(test_kerp,'check_') )
	suites.append( unittest.makeSuite(test_kei_zeros,'check_') )
	suites.append( unittest.makeSuite(test_keip_zeros,'check_') )
#F	suites.append( unittest.makeSuite(test_kelvin_zeros,'check_') )
	suites.append( unittest.makeSuite(test_ker_zeros,'check_') )
#E	suites.append( unittest.makeSuite(test_kerp_zeros,'check_') )	     
	suites.append( unittest.makeSuite(test_kn,'check_') )
	suites.append( unittest.makeSuite(test_kv,'check_') )
	suites.append( unittest.makeSuite(test_kve,'check_') )
#F	suites.append( unittest.makeSuite(test_kvp,'check_') )
	suites.append( unittest.makeSuite(test_laguerre,'check_') )
	suites.append( unittest.makeSuite(test_legendre,'check_') )
	suites.append( unittest.makeSuite(test_lmbda,'check_') )
	suites.append( unittest.makeSuite(test_log1p,'check_') )	
	suites.append( unittest.makeSuite(test_lpmn,'check_') )
#E	suites.append( unittest.makeSuite(test_lpn,'check_') )
	suites.append( unittest.makeSuite(test_lpmv,'check_') )
#!	suites.append( unittest.makeSuite(test_lqmn,'check_') )
	suites.append( unittest.makeSuite(test_lqn,'check_') )
###	 suites.append( unittest.makeSuite(test_mathieu_a,'check_') )
###	 suites.append( unittest.makeSuite(test_mathieu_even_coef,'check_') )
###	 suites.append( unittest.makeSuite(test_mathieu_odd_coef,'check_') )
###	 suites.append( unittest.makeSuite(test_modfresnelp,'check_') )
###	 suites.append( unittest.makeSuite(test_modfresnelm,'check_') )
#F	suites.append( unittest.makeSuite(test_obl_cv_seq,'check_') )	     
	suites.append( unittest.makeSuite(test_pbdn_seq,'check_') )
	suites.append( unittest.makeSuite(test_pbdv,'check_') )
#!	suites.append( unittest.makeSuite(test_pbdv_seq,'check_') )
###	suites.append( unittest.makeSuite(test_pbvv_seq,'check_') )
#E	suites.append( unittest.makeSuite(test_polygamma,'check_') )
#F	suites.append( unittest.makeSuite(test_pro_cv_seq,'check_') )
	suites.append( unittest.makeSuite(test_psi,'check_') )
	suites.append( unittest.makeSuite(test_radian,'check_') )	   
#F	suites.append( unittest.makeSuite(test_reshape,'check_') )
	suites.append( unittest.makeSuite(test_rgamma,'check_') )
	suites.append( unittest.makeSuite(test_riccati_jn,'check_') )
#F	suites.append( unittest.makeSuite(test_riccati_yn,'check_') )
#F	suites.append( unittest.makeSuite(test_round,'check_') )
	suites.append( unittest.makeSuite(test_sh_legendre,'check_') )
#F	suites.append( unittest.makeSuite(test_sh_chebyt,'check_') )
#F	suites.append( unittest.makeSuite(test_sh_chebyu,'check_') )
#F	suites.append( unittest.makeSuite(test_sh_jacobi,'check_') )	   
	suites.append( unittest.makeSuite(test_sinc,'check_') )
	suites.append( unittest.makeSuite(test_sindg,'check_') )
###	suites.append( unittest.makeSuite(test_sph_harm,'check_') )
#F	suites.append( unittest.makeSuite(test_sph_in,'check_') )
#E	suites.append( unittest.makeSuite(test_sph_inkn,'check_') )
#F	suites.append( unittest.makeSuite(test_sph_jn,'check_') )
#E	suites.append( unittest.makeSuite(test_sph_jnyn,'check_') )
#F	suites.append( unittest.makeSuite(test_sph_kn,'check_') )	   
	suites.append( unittest.makeSuite(test_sph_yn,'check_') )
	suites.append( unittest.makeSuite(test_take,'check_') )
	suites.append( unittest.makeSuite(test_tandg,'check_') )
	suites.append( unittest.makeSuite(test_y0,'check_') )
	suites.append( unittest.makeSuite(test_y1,'check_') )
#E	suites.append( unittest.makeSuite(test_y0_zeros,'check_') )
#F	suites.append( unittest.makeSuite(test_y1_zeros,'check_') )
#F	suites.append( unittest.makeSuite(test_y1p_zeros,'check_') )
	suites.append( unittest.makeSuite(test_yn_zeros,'check_') )
#F	suites.append( unittest.makeSuite(test_ynp_zeros,'check_') )	   
	suites.append( unittest.makeSuite(test_yn,'check_') )
	suites.append( unittest.makeSuite(test_yv,'check_') )
	suites.append( unittest.makeSuite(test_yve,'check_') )
	suites.append( unittest.makeSuite(test_yvp,'check_') )
	suites.append( unittest.makeSuite(test_zeros,'check_') )
        suites.append( unittest.makeSuite(test_general_function,'check_') )
    if level > 5:
	pass
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10,verbosity=1):
    all_tests = test_suite(level=level)
    runner = unittest.TextTestRunner(verbosity=verbosity)
    runner.run(all_tests)
    return runner

if __name__ == "__main__":
    test(1,2)
