"""
    This module is used for spline interpolation, and functions
    as a wrapper around the FITPACK Fortran interpolation
    package.  Its functionality is contained in the Spline
    class.

    Spline is primarily meant to be called by Interpolate1d
    or interp1d, but is a stand-alone class in its own right
    that is not available through these interfaces.
    
    Spline2d is primarily meant to be called by Interpolate2d
    or interp2d, but is a stand-alone class in its own right
    that is not available through these interfaces.

    The code has been modified from an older version of
    scipy.interpolate, where it was directly called by the
    user.  As such, it includes functionality not available through
    Interpolate1d.  For this reason, users may wish to get
    under the hood.

"""

import numpy as np
import warnings

import _dfitpack # extension module containing FITPACK subroutines in Fortran


class Spline(object):
    """ Univariate spline s(x) of degree k on the interval
        [xb,xe] calculated from a given set of data points
        (x,y).

        Can include least-squares fitting.

    """

    def __init__(self, x=None, y=None, w=None, bbox = [None]*2, k=3, s=0.0):
        """
        Input:
            x,y   - 1-d sequences of data points (x must be
                  in strictly ascending order)

        Optional input:
            k=3     - degree of the univariate spline.
            w       - positive 1-d sequence of weights
            bbox   - 2-sequence specifying the boundary of
                       the approximation interval.
                       By default, bbox=[x[0],x[-1]]
            s        - positive smoothing factor defined for
                       estimation condition:
                         sum((w[i]*( y[i]-s(x[i]) ))**2,axis=0) <= s
                        Default s=0
        """
        
        self._k = k
        self._s = s
        self._bbox = bbox
        self._w = w
        
        if x is not None and y is not None:
            self.init_xy(x, y)
            self._is_initialized = True
        else:
            self._is_initialized = False
        
    def init_xy(self, x, y):
        
        #_data == x,y,w,xb,xe,k,s,n,t,c,fp,fpint,nrdata,ier
        data = _dfitpack.fpcurf0(x, y, self._k, w=self._w,
                                xb=self._bbox[0], xe=self._bbox[1], s=self._s)
        if data[-1]==1:
            # nest too small, setting to maximum bound
            data = self._reset_nest(data)
        self._data = data
        # the relevant part of self._reset_class()
        n,t,c,k,ier = data[7],data[8],data[9],data[5],data[-1]
        self._eval_args = t[:n],c[:n],k
        
        self._is_initialized = True

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
        data = _dfitpack.fpcurf1(*args)
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
        data = _dfitpack.fpcurf1(*args)
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
        if self._is_initialized:
            if len(x) == 0: return np.array([]) #hack to cope with shape (0,)
            if nu is None:
                return _dfitpack.splev(*(self._eval_args+(x,)))
            return _dfitpack.splder(nu=nu,*(self._eval_args+(x,)))
        else:
            raise TypeError, "x and y must be set before interpolation is possible"

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
        return _dfitpack.splint(*(self._eval_args+(a,b)))

    def derivatives(self, x):
        """ Return all derivatives of the spline at the point x."""
        d,ier = _dfitpack.spalde(*(self._eval_args+(x,)))
        assert ier==0,`ier`
        return d

    def roots(self):
        """ Return the zeros of the spline.

        Restriction: only cubic splines are supported by fitpack.
        """
        k = self._data[5]
        if k==3:
            z,m,ier = _dfitpack.sproot(*self._eval_args[:2])
            assert ier==0,`ier`
            return z[:m]
        raise NotImplementedError,\
              'finding roots unsupported for non-cubic splines'
              

############################
## BELOW THIS POINT IS CODE FOR 2D INTERPOLATION
############################

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

        More commenting needed
        
        If (xi, yi) is outside the interpolation range, it is
        assigned the value of the nearest point that is within the
        interpolation range.
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
        """ Evaluate spline at positions (x[i], y[i]).
            x and y should be 1d arrays.
            
            If (xi, yi) is outside the interpolation range, it will be
            assigned the value of the nearest point that is within the
            interpolation range.
        """
        # FIXME : this function calls self.get_grid, which is extremely inefficient.  However,
        # I don't believe Fitpack really provides functionality to interpolate at scattered values.
        
        if self._is_initialized is False:
            raise Error, "x, y and z must be initialized before interpolating"
            
        # check input format
        assert ( isinstance(x, np.ndarray) and isinstance(y, np.ndarray) ), \
                    "newx and newy must both be numpy arrays"
        assert len(x) == len(y), "newx and newy must be of the same length"
        
        # sort only once for efficiency
        sorted_x = x.copy()
        sorted_x.sort()
        sorted_y = y.copy()
        sorted_y.sort()
        
        data_grid = self.get_grid(sorted_x, sorted_y)
        
        # FIXME : no list comprehension
        z = np.array([ data_grid[np.searchsorted(sorted_x, x[i]), np.searchsorted(sorted_y,y[i])] \
                                    for i,xi in enumerate(x) ])
        
        return z
        
        
    def get_grid(self, x, y):
        """ Evaluate spline at positions x[i],y[j]."""
        
        if self._is_initialized is not True:
            raise Error, "x, y and z must be initialized before interpolating"
        
        # check input format
        assert isinstance(x, np.ndarray) and isinstance(y, np.ndarray), \
                    "newx and newy must both be numpy arrays"
        assert len(x) == len(y), "newx and newy must be of the same length"
        
        tx,ty,c = self.tck[:3]
        kx,ky = self.degrees
        z,ier = _dfitpack.bispev(tx,ty,c,kx,ky,x,y)
        assert ier==0,'Invalid input: ier='+`ier`
        return z
        
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