import warnings
from numpy import zeros, concatenate, alltrue, ravel, all, diff
import numpy as np

import _dfitpack

_surfit_messages = {1:"""
        The required storage space exceeds the available storage space: nxest
        or nyest too small, or s too small.
        The weighted least-squares spline corresponds to the current set of
        knots.""",
                            2:"""
        A theoretically impossible result was found during the iteration
        process for finding a smoothing spline with fp = s: s too small or
        badly chosen eps.
        Weighted sum of squared residuals does not satisfy abs(fp-s)/s < tol.""",
                            3:"""
        the maximal number of iterations maxit (set to 20 by the program)
        allowed for finding a smoothing spline with fp=s has been reached:
        s too small.
        Weighted sum of squared residuals does not satisfy abs(fp-s)/s < tol.""",
                            4:"""
        No more knots can be added because the number of b-spline coefficients
        (nx-kx-1)*(ny-ky-1) already exceeds the number of data points m:
        either s or m too small.
        The weighted least-squares spline corresponds to the current set of
        knots.""",
                            5:"""
        No more knots can be added because the additional knot would (quasi)
        coincide with an old one: s too small or too large a weight to an
        inaccurate data point.
        The weighted least-squares spline corresponds to the current set of
        knots.""",
                            10:"""
        Error on entry, no approximation returned. The following conditions
        must hold:
        xb<=x[i]<=xe, yb<=y[i]<=ye, w[i]>0, i=0..m-1
        If iopt==-1, then
          xb<tx[kx+1]<tx[kx+2]<...<tx[nx-kx-2]<xe
          yb<ty[ky+1]<ty[ky+2]<...<ty[ny-ky-2]<ye""",
                            -3:"""
        The coefficients of the spline returned have been computed as the
        minimal norm least-squares solution of a (numerically) rank deficient
        system (deficiency=%i). If deficiency is large, the results may be
        inaccurate. Deficiency may strongly depend on the value of eps."""
                    }

class Spline2d(object):
    """ Bivariate spline s(x,y) of degrees kx and ky on the rectangle
        [xb,xe] x [yb, ye] calculated from a given set of data points
        (x,y,z).

        See also:

        bisplrep, bisplev - an older wrapping of FITPACK
        UnivariateSpline - a similar class for univariate spline interpolation
        SmoothUnivariateSpline - to create a BivariateSpline through the
                                 given points
        LSQUnivariateSpline - to create a BivariateSpline using weighted
                              least-squares fitting
    """
    def __init__(self, x=None, y=None, z=None, w=None, bbox=[None]*4, kx=3, ky=3, s=0.0, eps=None):
        """
            Input:
              x,y,z  - 1-d sequences of data points (order is not
                       important)
            Optional input:
              w          - positive 1-d sequence of weights
              bbox       - 4-sequence specifying the boundary of
                           the rectangular approximation domain.
                           By default, bbox=[min(x,tx),max(x,tx),
                                             min(y,ty),max(y,ty)]
              kx,ky=3,3  - degrees of the bivariate spline.
              s          - positive smoothing factor defined for
                           estimation condition:
                             sum((w[i]*(z[i]-s(x[i],y[i])))**2,axis=0) <= s
                           Default s=len(w) which should be a good value
                           if 1/w[i] is an estimate of the standard
                           deviation of z[i].
              eps        - a threshold for determining the effective rank
                           of an over-determined linear system of
                           equations. 0 < eps < 1, default is 1e-16.
        """
        
        self._w = w
        self._bbox = bbox
        self._kx = kx
        self._ky = kx
        self._s = s
        self._eps = eps
        
        if x is not None and y is not None and z is not None:
            self.init_xyz(x, y, z)
            self._is_initialized = True
        else:
            self._is_initialized = False
        
    def init_xyz(self, x, y, z):
        xb,xe,yb,ye = self._bbox
        nx,tx,ny,ty,c,fp,wrk1,ier = _dfitpack.surfit_smth(x,y,z,
                                                         self._w,
                                                         xb, xe, yb, ye,
                                                         self._kx, self._ky,
                                                         s=self._s,
                                                         eps=self._eps, lwrk2=1)
        if ier in [0,-1,-2]: # normal return
            pass
        else:
            message = _surfit_messages.get(ier,'ier=%s' % (ier))
            warnings.warn(message)

        self.fp = fp
        self.tck = tx[:nx],ty[:ny],c[:(nx-self._kx-1)*(ny-self._ky-1)]
        self.degrees = self._kx,self._ky
        
        self._is_initialized = True
        
    def __call__(self, x, y):
        """ Evaluate spline at positions x[i],y[i].
            x and y should be 1d arrays.
        """
        # what happens when x contains duplicate values?
        
        if self._is_initialized is not True:
            raise Error, "x, y and z must be initialized before interpolating"
        
        # sort only once for efficiency
        sorted_x = sorted(x)
        sorted_y = sorted(y)
        
        data_grid = self.get_grid(sorted_x, sorted_y)
        
        # fixme : no list comprehension
        z = np.array([ data_grid[np.searchsorted(sorted(x), x[i]), np.searchsorted(sorted(y),y[i])] \
                                    for i,xi in enumerate(x) ])
            
        return z
        
        
    def get_grid(self, x, y, mth='array'):
        """ Evaluate spline at positions x[i],y[j]."""
        
        if self._is_initialized is not True:
            raise Error, "x, y and z must be initialized before interpolating"
        
        if mth=='array':
            tx,ty,c = self.tck[:3]
            kx,ky = self.degrees
            z,ier = _dfitpack.bispev(tx,ty,c,kx,ky,x,y)
            assert ier==0,'Invalid input: ier='+`ier`
            return z
        raise NotImplementedError

    def get_residual(self):
        """ Return weighted sum of squared residuals of the spline
        approximation: sum ((w[i]*(z[i]-s(x[i],y[i])))**2,axis=0)
        """
        return self.fp
    def get_knots(self):
        """ Return a tuple (tx,ty) where tx,ty contain knots positions
            of the spline with respect to x-, y-variable, respectively.
            The position of interior and additional knots are given as
              t[k+1:-k-1] and t[:k+1]=b, t[-k-1:]=e, respectively.
        """
        return self.tck[:2]
    def get_coeffs(self):
        """ Return spline coefficients."""
        return self.tck[2]
    
    
    def integral(self, xa, xb, ya, yb):
        """
            Evaluate the integral of the spline over area [xa,xb] x [ya,yb].
            
            Parameters
            ----------
            xa, xb : float
                The end-points of the x integration interval.
            ya, yb : float
                The end-points of the y integration interval.
            
            Returns
            -------
            integ : float
                The value of the resulting integral.
            
        """
        tx,ty,c = self.tck[:3]
        kx,ky = self.degrees
        return _dfitpack.dblint(tx,ty,c,kx,ky,xa,xb,ya,yb)
        
# RectBivariateSpline in scipy.interpolate is for a rectangular grid and presumably must faster.  There are 3 levels of niceness: scattered
#       data, irregular grid, and regular grids.  Spline2d is for the first level and thus slow.  ndimage is for the 3rd level and thus fast.
#       I vote to no explicitly treat the 3rd level, but RecBivariateSpline does that if we want to implement it in the future.