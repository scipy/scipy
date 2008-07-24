"""
This module is used for spline interpolation, and functions
as a wrapper around the FITPACK Fortran interpolation
package.  It is not intended to be directly accessed by
the user, but rather through the class Interpolate1d.

The code has been modified from an older version of
scipy.interpolate, where it was directly called by the
user.  As such, it includes functionality not available through
Interpolate1d.  For this reason, users may wish to get
under the hood.

"""

import numpy as np

import dfitpack # lower-level wrapper around FITPACK


class Spline(object):
    """ Univariate spline s(x) of degree k on the interval
    [xb,xe] calculated from a given set of data points
    (x,y).

    Can include least-squares fitting.

    See also:

    splrep, splev, sproot, spint, spalde - an older wrapping of FITPACK
    BivariateSpline - a similar class for bivariate spline interpolation
    """

    def __init__(self, x, y, w=None, bbox = [None]*2, k=3, s=0.0):
        """
        Input:
          x,y   - 1-d sequences of data points (x must be
                  in strictly ascending order)

        Optional input:
          w          - positive 1-d sequence of weights
          bbox       - 2-sequence specifying the boundary of
                       the approximation interval.
                       By default, bbox=[x[0],x[-1]]
          k=3        - degree of the univariate spline.
          s          - positive smoothing factor defined for
                       estimation condition:
                         sum((w[i]*( y[i]-s(x[i]) ))**2,axis=0) <= s
                       Default s=len(w) which should be a good value
                       if 1/w[i] is an estimate of the standard
                       deviation of y[i].
        """
        #_data == x,y,w,xb,xe,k,s,n,t,c,fp,fpint,nrdata,ier
        data = dfitpack.fpcurf0(x, y, k, w=w,
                                xb=bbox[0], xe=bbox[1], s=s)
        if data[-1]==1:
            # nest too small, setting to maximum bound
            data = self._reset_nest(data)
        self._data = data
        # the relevant part of self._reset_class()
        n,t,c,k,ier = data[7],data[8],data[9],data[5],data[-1]
        self._eval_args = t[:n],c[:n],k

    def _reset_nest(self, data, nest=None):
        n = data[10]
        if nest is None:
            k,m = data[5],len(data[0])
            nest = m+k+1 # this is the maximum bound for nest
        else:
            assert n<=nest,"nest can only be increased"
        t,c,fpint,nrdata = data[8].copy(),data[9].copy(),\
                           data[11].copy(),data[12].copy()
        t.resize(nest)
        c.resize(nest)
        fpint.resize(nest)
        nrdata.resize(nest)
        args = data[:8] + (t,c,n,fpint,nrdata,data[13])
        data = dfitpack.fpcurf1(*args)
        return data

    def set_smoothing_factor(self, s):
        """ Continue spline computation with the given smoothing
        factor s and with the knots found at the last call.
        
        """
        data = self._data
        if data[6]==-1:
            warnings.warn('smoothing factor unchanged for'
                          'LSQ spline with fixed knots')
            return
        args = data[:6] + (s,) + data[7:]
        data = dfitpack.fpcurf1(*args)
        if data[-1]==1:
            # nest too small, setting to maximum bound
            data = self._reset_nest(data)
        self._data = data
        #self._reset_class()

    def __call__(self, x, nu=None):
        """ Evaluate spline (or its nu-th derivative) at positions x.
        Note: x can be unordered but the evaluation is more efficient
        if x is (partially) ordered.
        
        """
        
        if len(x) == 0: return np.array([]) #hack to cope with shape (0,)
        if nu is None:
            return dfitpack.splev(*(self._eval_args+(x,)))
        return dfitpack.splder(nu=nu,*(self._eval_args+(x,)))

    def get_knots(self):
        """ Return the positions of (boundary and interior)
        knots of the spline.
        """
        data = self._data
        k,n = data[5],data[7]
        return data[8][k:n-k]

    def get_coeffs(self):
        """Return spline coefficients."""
        data = self._data
        k,n = data[5],data[7]
        return data[9][:n-k-1]

    def get_residual(self):
        """Return weighted sum of squared residuals of the spline
        approximation: sum ((w[i]*(y[i]-s(x[i])))**2,axis=0)
        
        """
        return self._data[10]

    def integral(self, a, b):
        """ Return definite integral of the spline between two
        given points.
        """
        return dfitpack.splint(*(self._eval_args+(a,b)))

    def derivatives(self, x):
        """ Return all derivatives of the spline at the point x."""
        d,ier = dfitpack.spalde(*(self._eval_args+(x,)))
        assert ier==0,`ier`
        return d

    def roots(self):
        """ Return the zeros of the spline.

        Restriction: only cubic splines are supported by fitpack.
        """
        k = self._data[5]
        if k==3:
            z,m,ier = dfitpack.sproot(*self._eval_args[:2])
            assert ier==0,`ier`
            return z[:m]
        raise NotImplementedError,\
              'finding roots unsupported for non-cubic splines'
                    


# testing
import unittest
import time
from numpy import arange, allclose, ones

class Test(unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(x, y))
        
    def test_linearSpl(self):
        N = 3000.
        x = np.arange(N)
        y = np.arange(N)
        T1 = time.clock()
        interp_func = Spline(x, y, k=1)
        T2 = time.clock()
        print "time to create order 1 spline interpolation function with N = %i:" % N, T2 - T1
        new_x = np.arange(N)+0.5
        t1 = time.clock()
        new_y = interp_func(new_x)
        t2 = time.clock()
        print "time for order 1 spline interpolation with N = %i:" % N, t2 - t1
        self.assertAllclose(new_y[:5], [0.5, 1.5, 2.5, 3.5, 4.5])
    
    def runTest(self):
        self.test_linearSpl()
                             
if __name__ == '__main__':
    unittest.main()
    
    