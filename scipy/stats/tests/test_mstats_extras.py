# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for maskedArray statistics.

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"

import numpy as np

import numpy.ma as ma

import scipy.stats.mstats as ms
#import scipy.stats.mmorestats as mms

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_almost_equal, assert_


class TestMisc(TestCase):
    #
    def __init__(self, *args, **kwargs):
        TestCase.__init__(self, *args, **kwargs)
    #
    def test_mjci(self):
        "Tests the Marits-Jarrett estimator"
        data = ma.array([ 77, 87, 88,114,151,210,219,246,253,262,
                          296,299,306,376,428,515,666,1310,2611])
        assert_almost_equal(ms.mjci(data),[55.76819,45.84028,198.87875],5)
    #
    def test_trimmedmeanci(self):
        "Tests the confidence intervals of the trimmed mean."
        data = ma.array([545,555,558,572,575,576,578,580,
                         594,605,635,651,653,661,666])
        assert_almost_equal(ms.trimmed_mean(data,0.2), 596.2, 1)
        assert_equal(np.round(ms.trimmed_mean_ci(data,(0.2,0.2)),1),
                     [561.8, 630.6])
    #
    def test_idealfourths(self):
        "Tests ideal-fourths"
        test = np.arange(100)
        assert_almost_equal(np.asarray(ms.idealfourths(test)),
                            [24.416667,74.583333],6)
        test_2D = test.repeat(3).reshape(-1,3)
        assert_almost_equal(ms.idealfourths(test_2D, axis=0),
                            [[24.416667,24.416667,24.416667],
                             [74.583333,74.583333,74.583333]],6)
        assert_almost_equal(ms.idealfourths(test_2D, axis=1),
                            test.repeat(2).reshape(-1,2))
        test = [0,0]
        _result = ms.idealfourths(test)
        assert_(np.isnan(_result).all())

#..............................................................................
class TestQuantiles(TestCase):
    #
    def __init__(self, *args, **kwargs):
        TestCase.__init__(self, *args, **kwargs)
    #
    def test_hdquantiles(self):
        data = [0.706560797,0.727229578,0.990399276,0.927065621,0.158953014,
            0.887764025,0.239407086,0.349638551,0.972791145,0.149789972,
            0.936947700,0.132359948,0.046041972,0.641675031,0.945530547,
            0.224218684,0.771450991,0.820257774,0.336458052,0.589113496,
            0.509736129,0.696838829,0.491323573,0.622767425,0.775189248,
            0.641461450,0.118455200,0.773029450,0.319280007,0.752229111,
            0.047841438,0.466295911,0.583850781,0.840581845,0.550086491,
            0.466470062,0.504765074,0.226855960,0.362641207,0.891620942,
            0.127898691,0.490094097,0.044882048,0.041441695,0.317976349,
            0.504135618,0.567353033,0.434617473,0.636243375,0.231803616,
            0.230154113,0.160011327,0.819464108,0.854706985,0.438809221,
            0.487427267,0.786907310,0.408367937,0.405534192,0.250444460,
            0.995309248,0.144389588,0.739947527,0.953543606,0.680051621,
            0.388382017,0.863530727,0.006514031,0.118007779,0.924024803,
            0.384236354,0.893687694,0.626534881,0.473051932,0.750134705,
            0.241843555,0.432947602,0.689538104,0.136934797,0.150206859,
            0.474335206,0.907775349,0.525869295,0.189184225,0.854284286,
            0.831089744,0.251637345,0.587038213,0.254475554,0.237781276,
            0.827928620,0.480283781,0.594514455,0.213641488,0.024194386,
            0.536668589,0.699497811,0.892804071,0.093835427,0.731107772]
        #
        assert_almost_equal(ms.hdquantiles(data,[0., 1.]),
                            [0.006514031, 0.995309248])
        hdq = ms.hdquantiles(data,[0.25, 0.5, 0.75])
        assert_almost_equal(hdq, [0.253210762, 0.512847491, 0.762232442,])
        hdq = ms.hdquantiles_sd(data,[0.25, 0.5, 0.75])
        assert_almost_equal(hdq, [0.03786954, 0.03805389, 0.03800152,], 4)
        #
        data = np.array(data).reshape(10,10)
        hdq = ms.hdquantiles(data,[0.25,0.5,0.75],axis=0)
        assert_almost_equal(hdq[:,0], ms.hdquantiles(data[:,0],[0.25,0.5,0.75]))
        assert_almost_equal(hdq[:,-1], ms.hdquantiles(data[:,-1],[0.25,0.5,0.75]))
        hdq = ms.hdquantiles(data,[0.25,0.5,0.75],axis=0,var=True)
        assert_almost_equal(hdq[...,0],
                            ms.hdquantiles(data[:,0],[0.25,0.5,0.75],var=True))
        assert_almost_equal(hdq[...,-1],
                            ms.hdquantiles(data[:,-1],[0.25,0.5,0.75], var=True))


###############################################################################

if __name__ == "__main__":
    run_module_suite()
