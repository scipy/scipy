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
    def _check_bdtrik(self):
        assert_equal(cephes.bdtrik(1,3,0.5),3.0)


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
    def _check_chndtrix(self):
        assert_equal(cephes.chndtrix(0,1,0),0.0)

    def check_gamma(self):
        assert_equal(cephes.gamma(5),24)
    def _check_ncfdtr(self):
        assert_equal(cephes.ncfdtr(2,2,0,1),0.5)
        
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
