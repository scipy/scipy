"""
fitpack --- curve and surface fitting with splines

fitpack is based on a collection of Fortran routines DIERCKX
by P. Dierckx (see http://www.netlib.org/dierckx/) transformed
to double routines by Pearu Peterson.
"""
# Created by Pearu Peterson, June,August 2003

__all__ = [
    'UnivariateSpline',
    'InterpolatedUnivariateSpline',
    'LSQUnivariateSpline',

    'LSQBivariateSpline',
    'SmoothBivariateSpline']

import warnings
from scipy_base import zeros, Float, concatenate, alltrue
from scipy_distutils.misc_util import PostponedException

try: import dfitpack
except: dfitpack = PostponedException()

################ Univariate spline ####################

_curfit_messages = {1:"""
The required storage space exceeds the available storage space, as
specified by the parameter nest: nest too small. If nest is already
large (say nest > m/2), it may also indicate that s is too small.
The approximation returned is the weighted least-squares spline
according to the knots t[0],t[1],...,t[n-1]. (n=nest) the parameter fp
gives the corresponding weighted sum of squared residuals (fp>s).
""",
                    2:"""
A theoretically impossible result was found during the iteration
proces for finding a smoothing spline with fp = s: s too small.
There is an approximation returned but the corresponding weighted sum
of squared residuals does not satisfy the condition abs(fp-s)/s < tol.""",
                    3:"""
The maximal number of iterations maxit (set to 20 by the program)
allowed for finding a smoothing spline with fp=s has been reached: s
too small.
There is an approximation returned but the corresponding weighted sum
of squared residuals does not satisfy the condition abs(fp-s)/s < tol.""",
                    10:"""
Error on entry, no approximation returned. The following conditions
must hold:
xb<=x[0]<x[1]<...<x[m-1]<=xe, w[i]>0, i=0..m-1
if iopt=-1:
  xb<t[k+1]<t[k+2]<...<t[n-k-2]<xe"""
                    }

class UnivariateSpline:
    """ Univariate spline s(x) of degree k on the interval
    [xb,xe] calculated from a given set of data points
    (x,y).
    """

    def __init__(self, x, y, w=None, bbox = [None]*2, k=3, s=None):
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
                         sum((w[i]*(y[i]-s(x[i])))**2) <= s
                       Default s=len(w) which should be a good value
                       if 1/w[i] is an estimate of the standard
                       deviation of y[i].
        """
        #_data == x,y,w,xb,xe,k,s,n,t,c,fp,fpint,nrdata,ier
        data = dfitpack.fpcurf0(x,y,k,w=w,
                                xb=bbox[0],xe=bbox[1],s=s)
        if data[-1]==1:
            # nest too small, setting to maximum bound
            data = self._reset_nest(data)
        self._data = data
        self._reset_class()

    def _reset_class(self):
        data = self._data
        n,t,c,k,ier = data[7],data[8],data[9],data[5],data[-1]
        self._eval_args = t[:n],c[:n-k-1],k
        if ier==0:
            # the spline returned has a residual sum of squares fp
            # such that abs(fp-s)/s <= tol with tol a relative
            # tolerance set to 0.001 by the program
            pass
        elif ier==-1:
            # the spline returned is an interpolating spline
            self.__class__ = InterpolatedUnivariateSpline
        elif ier==-2:
            # the spline returned is the weighted least-squares
            # polynomial of degree k. In this extreme case fp gives
            # the upper bound fp0 for the smoothing factor s.
            self.__class__ = LSQUnivariateSpline
        else:
            # error
            if ier==1:
                self.__class__ = LSQUnivariateSpline
            message = _curfit_messages.get(ier,'ier=%s' % (ier))
            warnings.warn(message)

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
        self._reset_class()

    def __call__(self, x, nu=None):
        """ Evaluate spline (or its nu-th derivative) at positions x.
        Note: x can be unordered but the evaluation is
        more efficient if x is (partially) ordered.
        """
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
        """ Return spline coefficients."""
        data = self._data
        k,n = data[5],data[7]
        return data[9][:n-k-1]

    def get_residual(self):
        """ Return weighted sum of squared residuals of the spline
        approximation: sum ((w[i]*(y[i]-s(x[i])))**2)
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
        raise NotImplementedError,'finding roots unsupported for non-cubic splines'

class InterpolatedUnivariateSpline(UnivariateSpline):
    """ Interpolated univariate spline approximation."""
    def __init__(self, x, y, w=None, bbox = [None]*2, k=3):
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
        """
        #_data == x,y,w,xb,xe,k,s,n,t,c,fp,fpint,nrdata,ier
        self._data = dfitpack.fpcurf0(x,y,k,w=w,
                                      xb=bbox[0],xe=bbox[1],s=0)
        self._reset_class()

class LSQUnivariateSpline(UnivariateSpline):
    """ Weighted least-squares univariate spline approximation."""

    def __init__(self, x, y, t, w=None, bbox = [None]*2, k=3):
        """
        Input:
          x,y   - 1-d sequences of data points (x must be
                  in strictly ascending order)
          t     - 1-d sequence of the positions of user-defined
                  interior knots of the spline (t must be in strictly
                  ascending order and bbox[0]<t[0]<...<t[-1]<bbox[-1])

        Optional input:
          w          - positive 1-d sequence of weights
          bbox       - 2-sequence specifying the boundary of
                       the approximation interval.
                       By default, bbox=[x[0],x[-1]]
          k=3        - degree of the univariate spline.
        """
        #_data == x,y,w,xb,xe,k,s,n,t,c,fp,fpint,nrdata,ier
        xb=bbox[0]
        xe=bbox[1]
        if xb is None: xb = x[0]
        if xe is None: xe = x[-1]
        t = concatenate(([xb]*(k+1),t,[xe]*(k+1)))
        n = len(t)
        if not alltrue(t[k+1:n-k]-t[k:n-k-1] > 0):
            raise ValueError,\
                  'Interior knots t must satisfy Schoenberg-Whitney conditions'
        data = dfitpack.fpcurfm1(x,y,k,t,w=w,xb=xb,xe=xe)
        self._data = data[:-3] + (None,None,data[-1])
        self._reset_class()


################ Bivariate spline ####################

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


class BivariateSpline:
    """ Bivariate spline s(x,y) of degrees kx and ky on the rectangle
    [xb,xe] x [yb, ye] calculated from a given set of data points
    (x,y,z).
    """

    def get_residual(self):
        """ Return weighted sum of squared residuals of the spline
        approximation: sum ((w[i]*(z[i]-s(x[i],y[i])))**2)
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
    def __call__(self,x,y,mth='array'):
        """ Evaluate spline at positions x,y."""
        if mth=='array':
            tx,ty,c = self.tck[:3]
            kx,ky = self.degrees
            z,ier = dfitpack.bispev(tx,ty,c,kx,ky,x,y)
            assert ier==0,'Invalid input: ier='+`ier`
            return z
        raise NotImplementedError

class SmoothBivariateSpline(BivariateSpline):
    """ Smooth bivariate spline approximation."""

    def __init__(self, x, y, z, w=None,
                 bbox = [None]*4, kx=3, ky=3, s=None, eps=None):
        """
        Input:
          x,y,z  - 1-d sequences of data points (order is not
                   important) 
        Optional input:
          w          - positive 1-d sequence of weights
          bbox       - 4-sequence specifying the boundary of
                       the rectangular approximation domain.
                       By default, bbox=[min(x,tx),max(x,tx),min(y,ty),max(y,ty)]
          kx,ky=3,3  - degrees of the bivariate spline.
          s          - positive smoothing factor defined for
                       estimation condition:
                         sum((w[i]*(z[i]-s(x[i],y[i])))**2) <= s
                       Default s=len(w) which should be a good value
                       if 1/w[i] is an estimate of the standard
                       deviation of z[i].
          eps        - a threshold for determining the effective rank
                       of an over-determined linear system of
                       equations. 0 < eps < 1, default is 1e-16.
        """
        xb,xe,yb,ye = bbox
        nx,tx,ny,ty,c,fp,wrk1,ier = dfitpack.surfit_smth(x,y,z,w,\
                                                         xb,xe,yb,ye,\
                                                         kx,ky,s=s,eps=eps,lwrk2=1)
        if ier in [0,-1,-2]: # normal return
            pass
        else:
            message = _surfit_messages.get(ier,'ier=%s' % (ier))
            warnings.warn(message)

        self.fp = fp
        self.tck = tx[:nx],ty[:ny],c[:(nx-kx-1)*(ny-ky-1)]
        self.degrees = kx,ky

class LSQBivariateSpline(BivariateSpline):
    """ Weighted least-squares spline approximation.
    """
    def __init__(self, x, y, z, tx, ty, w=None,
                 bbox = [None]*4,
                 kx=3, ky=3, eps=None):
        """
        Input:
          x,y,z  - 1-d sequences of data points (order is not
                   important) 
          tx,ty  - strictly ordered 1-d sequences of knots
                   coordinates. 
        Optional input:
          w          - positive 1-d sequence of weights
          bbox       - 4-sequence specifying the boundary of
                       the rectangular approximation domain.
                       By default, bbox=[min(x,tx),max(x,tx),min(y,ty),max(y,ty)]
          kx,ky=3,3  - degrees of the bivariate spline.
          eps        - a threshold for determining the effective rank
                       of an over-determined linear system of
                       equations. 0 < eps < 1, default is 1e-16. 
        """
        nx = 2*kx+2+len(tx)
        ny = 2*ky+2+len(ty)
        tx1 = zeros((nx,),Float)
        ty1 = zeros((ny,),Float)
        tx1[kx+1:nx-kx-1] = tx
        ty1[ky+1:ny-ky-1] = ty

        xb,xe,yb,ye = bbox
        tx1,ty1,c,fp,ier = dfitpack.surfit_lsq(x,y,z,tx1,ty1,w,\
                                               xb,xe,yb,ye,\
                                               kx,ky,eps,lwrk2=1)
        if ier>10:
            tx1,ty1,c,fp,ier = dfitpack.surfit_lsq(x,y,z,tx1,ty1,w,\
                                                   xb,xe,yb,ye,\
                                                   kx,ky,eps,lwrk2=ier)
        if ier in [0,-1,-2]: # normal return
            pass
        else:
            if ier<-2:
                deficiency = (nx-kx-1)*(ny-ky-1)+ier
                message = _surfit_messages.get(-3) % (deficiency)
            else:
                message = _surfit_messages.get(ier,'ier=%s' % (ier))
            warnings.warn(message)
        self.fp = fp
        self.tck = tx1,ty1,c
        self.degrees = kx,ky

