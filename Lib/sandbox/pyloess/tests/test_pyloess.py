"""
Wrapper to lowess and stl routines.

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date: 2007-02-28 02:23:25 -0500 (Wed, 28 Feb 2007) $
:version: $Id: generic.py 145 2007-02-28 07:23:25Z backtopop $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 145 $"
__date__     = '$Date: 2007-02-28 02:23:25 -0500 (Wed, 28 Feb 2007) $'

import os

import numpy
from numpy import bool_, complex_, float_, int_, str_, object_
import numpy.core.numeric as numeric
from numpy.core.records import recarray

from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg, \
        assert_equal, assert_almost_equal

import pyloess
reload(pyloess)
from pyloess import lowess, stl


def get_co2data():
    "Reads CO2 data."
    filename = os.path.join('tests','co2_data')
    F = open(filename, 'r')
    data = []
    for line in F.readlines():
        data.append([float(x) for x in line.rstrip().split()])
    return numpy.concatenate(data)

def get_co2results():
    "Gets theoretical results of smoothed CO2."
    filename = os.path.join('tests','co2_results_double')
    F = open(filename, 'r')
    result = []
    for line in F.readlines():
        result.append(numpy.fromiter((float(x) for x in line.rstrip().split()),
                                     float_))
    return result

def set_parameters():
    "Returns the parameters of the STL on CO2 data."
    parameters = dict(np=12, ns=35, nt=19, nl=13, no=2, ni=1,
                      nsjump=4, ntjump=2, nljump=2,
                      isdeg=1, itdeg=1, ildeg=1)
    return parameters
########################################################################
class test_lowess(NumpyTestCase):
    "Test class for lowess."
    #
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        X = [ 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8,10,12,14,50]
        Y = [18, 2,15, 6,10, 4,16,11, 7, 3,14,17,20,12, 9,13, 1, 8, 5,19]
        self.data = (X, Y)
    #
    def test_lowess_1(self):
        "Tests lowess on typical data. part #1."
        (X, Y) = self.data
        YS = [13.659,11.145, 8.701, 9.722,10.000,11.300,11.300,11.300,
              11.300,11.300,11.300,11.300,11.300,11.300,11.300,13.000, 
              6.440, 5.596,  5.456,18.998]
        Z = lowess(X, Y, f=0.25, nsteps=0, delta=0)
        assert_almost_equal(Z.smooth, YS, decimal=3)
        assert_almost_equal(Z.residuals+Z.smooth, Y, decimal=3)
    #
    def test_lowess_2(self):
        "Tests lowess on typical data. part #2."
        (X, Y) = self.data
        YS = [13.659,12.347,11.034, 9.722,10.511,11.300,11.300,11.300,
              11.300,11.300,11.300,11.300,11.300,11.300,11.300,13.000, 
               6.440, 5.596, 5.456,18.998]
        Z = lowess(X, Y, f=0.25, nsteps=0, delta=3)
        assert_almost_equal(Z.smooth, YS, decimal=3)
        assert_almost_equal(Z.residuals+Z.smooth, Y, decimal=3)        
    #
    def test_lowess_3(self):
        "Tests lowess on typical data. part #3."
        (X, Y) = self.data
        YS = [14.811,12.115, 8.984, 9.676,10.000,11.346,11.346,11.346,
              11.346,11.346,11.346,11.346,11.346,11.346,11.346,13.000, 
               6.734, 5.744, 5.415,18.998 ]
        Z = lowess(X, Y, f=0.25, nsteps=2, delta=0)
        assert_almost_equal(Z.smooth, YS, decimal=3)
        assert_almost_equal(Z.residuals+Z.smooth, Y, decimal=3)        

#-----------------------------------------------------------------------
class test_stl(NumpyTestCase):
    "Tests STL."
    #
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        co2_data = get_co2data()
        co2_results = get_co2results()
        parameters = set_parameters()
        self.d = (co2_data, co2_results, parameters)
    #
    def test_stl_1(self):
        "Tests a classic STL."
        (co2_data, co2_results, parameters) = self.d
        co2_fitted = stl(co2_data, robust=False, **parameters)
        assert_almost_equal(co2_fitted.seasonal, co2_results[0], 6)
        assert_almost_equal(co2_fitted.trend, co2_results[1], 6)
        assert_almost_equal(co2_fitted.weights, co2_results[2], 6)
    #
    def test_stl_2(self):
        "Tests a robust STL."
        (co2_data, co2_results, parameters) = self.d
        co2_fitted = stl(co2_data, robust=True, **parameters)
        assert_almost_equal(co2_fitted.seasonal, co2_results[4], 6)
        assert_almost_equal(co2_fitted.trend, co2_results[5], 6)
        assert_almost_equal(co2_fitted.weights, co2_results[6], 6)
        
########################################################################
if __name__ == '__main__':
    NumpyTest().run()       
    #
    co2_data = get_co2data()
    co2_results = get_co2results()
    parameters = set_parameters()
    co2_fitted = stl(co2_data, robust=False, **parameters)
    

