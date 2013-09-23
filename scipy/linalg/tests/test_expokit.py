#!/usr/bin/env python 
#
# Created by: Alex Leach, November 2012
#

import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_array_almost_equal, \
    assert_array_almost_equal_nulp

from scipy.linalg.expokit import expm

#from scipy.linalg.expokit import dgchbv, dgexpv, dgpadm, dgphiv, dmexpv, dnchbv, \
#    dschbv, dsexpv, dspadm, dsphiv, zgchbv, zgexpv, zgpadm, zgphiv, zhexpv,  \
#    zhpadm, zhphiv, znchbv


class TestPADM(TestCase):
    def test_zero(self):
        a = array([[0.,0],[0,0]])
        assert_array_almost_equal(expm(a),[[1,0],[0,1]])

    def test_simple(self):
        a = array([[0.,1.],[0,0]])
        assert_array_almost_equal(expm(a),[[1,1],[0,1]] )

    def test_complex(self):
        a = array([[0.,1.],[0,0]])
        assert_array_almost_equal(expm(a),[[1,1],[0,1]] )

if __name__ == '__main__':
    run_module_suite()
