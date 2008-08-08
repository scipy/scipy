#!/usr/bin/env python
# Created by Pearu Peterson, June 2003
# Modified by John Travers, October 2006
""" Test functions for spline.spline module
"""
__usage__ = """
Build spline:
  python setup_spline.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.spline.test(<level>)'
Run tests if spline is not installed:
  python tests/test_spline.py [<level>]
"""

import sys
from numpy.testing import *
from numpy import array, arange, around, pi, sin, cos


from scipy.sandbox.spline.spline import UnivariateSpline,LSQUnivariateSpline,\
     InterpolatedUnivariateSpline, LSQBivariateSpline, SmoothBivariateSpline,\
     RectBivariateSpline

from dierckx_test_data import *

class TestUnivariateSpline(TestCase):
    def test_linear_constant(self):
        x = [1,2,3]
        y = [3,3,3]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[3,3,3])

    def test_linear_1d(self):
        x = [1,2,3]
        y = [0,2,4]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[0,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[0,1,2])

    def test_curfit_against_dierckx(self):
        """ Test against results obtined from the pure fortran routines.

            Here we check simple spline creation and evaluation.
        """
        x,y = curfit_test['x'],curfit_test['y']
        k,s = curfit_test_smth['k'],curfit_test_smth['s']
        iopt = curfit_test_smth['iopt']
        for i in range(len(k)):
            if iopt[i] == 0:
                uspl = UnivariateSpline(x,y,k=k[i],s=s[i])
            elif iopt[i] == 1:
                uspl.set_smoothing_factor(s[i])
            assert_almost_equal(uspl.get_residual(),
                                curfit_test_smth['fp'][i], decimal=2)
            n = uspl._data[7]
            assert_equal(n,curfit_test_smth['n'][i])
            assert_array_almost_equal(around(uspl.get_knots(),1),
                                      curfit_test_smth['t'][i][k[i]:n-k[i]])
            assert_array_almost_equal(around(uspl.get_coeffs(),4),
                                      curfit_test_smth['c'][i], decimal=3)
            assert_array_almost_equal(around(uspl(x),1),
                                      curfit_test_smth['sp'][i])

    def test_spint_spalde(self):
        per = [0, 0, 0]
        N = [20, 20, 50]
        ia = [0, 0, 0.2*pi]
        ib = [0, 0, pi]
        a,b = 0,2*pi
        dx = 0.2*pi
        k = range(1,6)
        for i in range(len(per)):
            x=a+(b-a)*arange(N[i]+1,dtype=float)/float(N[i])
            v=f1(x)
            for j in range(len(k)):
                uspl = UnivariateSpline(x,v,k=k[j],s=0)
                ir = uspl.integral(ia[i],ib[i])
                dr = uspl.derivatives(dx)
                assert_almost_equal(ir, f1(ib[i],-1)-f1(ia[i],-1), decimal=2)
                d=0
                for ddr in dr:
                    if d<k[j]-1:
                        assert_almost_equal(1, ddr/f1(dx,d), decimal=2)
                    d=d+1

    def test_sproot(self):
        a=0
        b=15
        N=20
        x=a+(b-a)*arange(N+1,dtype=float)/float(N)
        v=f1(x)
        k=3
        uspl = UnivariateSpline(x,v,k=k,s=0)
        ex = array([0.0, pi, 2.0*pi, 3.0*pi, 4.0*pi])
        assert_array_almost_equal(uspl.roots(),ex, decimal=3)

class TestLSQUnivariateSpline(TestCase):
    def test_curfit_against_dierckx(self):
        """ Test against results obtined from the pure fortran routines.

            Here we check simple spline creation and evaluation.
        """
        x,y = curfit_test['x'],curfit_test['y']
        k = curfit_test_lsq['k']
        for i in range(len(k)):
            t = curfit_test_lsq['t'][i]
            lsquspl = LSQUnivariateSpline(x,y,t,k=k[i])
            assert_almost_equal(lsquspl.get_residual(),
                                curfit_test_lsq['fp'][i], decimal=2)
            assert_array_almost_equal(around(lsquspl.get_coeffs(),4),
                                      curfit_test_lsq['c'][i], decimal=3)
            assert_array_almost_equal(around(lsquspl(x),1),
                                      curfit_test_lsq['sp'][i])

class TestLSQBivariateSpline(TestCase):
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)
        #print lut.get_knots()
        #print lut.get_coeffs()
        #print lut.get_residual()

class TestSmoothBivariateSpline(TestCase):
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[3,3,3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[3,3],[3,3],[3,3]])

    def test_linear_1d(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,0,0,2,2,2,4,4,4]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[0,0,4,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[0,0],[1,1],[2,2]])

class TestRectBivariateSpline(TestCase):
    def test_defaults(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y),z)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
