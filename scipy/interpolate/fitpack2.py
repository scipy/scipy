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

    'BivariateSpline',
    'LSQBivariateSpline',
    'SmoothBivariateSpline',
    'RectBivariateSpline']

import warnings
from numpy import zeros, concatenate, alltrue, ravel, all, diff

import dfitpack

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

class UnivariateSpline(object):
    """
    One-dimensional smoothing spline fit to a given set of data points.

    Fits a spline y=s(x) of degree `k` to the provided `x`,`y` data.  `s`
    specifies the number of knots by specifying a smoothing condition.

    Parameters
    ----------
    x : sequence
        input dimension of data points -- must be increasing
    y : sequence
         input dimension of data points
    w : sequence or None, optional
        weights for spline fitting.  Must be positive.  If None (default),
        weights are all equal.
    bbox : sequence or None, optional
        2-sequence specifying the boundary of the approximation interval. If
        None (default), bbox=[x[0],x[-1]].
    k : int, optional
        Degree of the smoothing spline.  Must be <= 5.
    s : float or None, optional
        Positive smoothing factor used to choose the number of knots.  Number
        of knots will be increased until the smoothing condition is satisfied:

        sum((w[i]*(y[i]-s(x[i])))**2,axis=0) <= s

        If None (default), s=len(w) which should be a good value if 1/w[i] is
        an estimate of the standard deviation of y[i].  If 0, spline will
        interpolate through all data points.


    See Also
    --------
    InterpolatedUnivariateSpline : Subclass with smoothing forced to 0
    LSQUnivariateSpline : Subclass in which knots are user-selected instead of
        being set by smoothing condition
    splrep : An older, non object-oriented wrapping of FITPACK
    splev, sproot, splint, spalde
    BivariateSpline : A similar class for two-dimensional spline interpolation



    Examples
    --------
    >>> from numpy import linspace,exp
    >>> from numpy.random import randn
    >>> from scipy.interpolate import UnivariateSpline
    >>> x = linspace(-3,3,100)
    >>> y = exp(-x**2) + randn(100)/10
    >>> s = UnivariateSpline(x,y,s=1)
    >>> xs = linspace(-3,3,1000)
    >>> ys = s(xs)

    xs,ys is now a smoothed, super-sampled version of the noisy gaussian x,y

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
                         sum((w[i]*(y[i]-s(x[i])))**2,axis=0) <= s
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
        self._eval_args = t[:n],c[:n],k
        if ier==0:
            # the spline returned has a residual sum of squares fp
            # such that abs(fp-s)/s <= tol with tol a relative
            # tolerance set to 0.001 by the program
            pass
        elif ier==-1:
            # the spline returned is an interpolating spline
            self._set_class(InterpolatedUnivariateSpline)
        elif ier==-2:
            # the spline returned is the weighted least-squares
            # polynomial of degree k. In this extreme case fp gives
            # the upper bound fp0 for the smoothing factor s.
            self._set_class(LSQUnivariateSpline)
        else:
            # error
            if ier==1:
                self._set_class(LSQUnivariateSpline)
            message = _curfit_messages.get(ier,'ier=%s' % (ier))
            warnings.warn(message)

    def _set_class(self, cls):
        self._spline_class = cls
        if self.__class__ in (UnivariateSpline, InterpolatedUnivariateSpline,
                              LSQUnivariateSpline):
            self.__class__ = cls
        else:
            # It's an unknown subclass -- don't change class. cf. #731
            pass

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
        Note: x can be unordered but the evaluation is more efficient
        if x is (partially) ordered.

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

class InterpolatedUnivariateSpline(UnivariateSpline):
    """
    One-dimensional interpolating spline for a given set of data points.

    Fits a spline y=s(x) of degree `k` to the provided `x`,`y` data. Spline
    function passes through all provided points. Equivalent to
    `UnivariateSpline` with  s=0.

    Parameters
    ----------
    x : sequence
        input dimension of data points -- must be increasing
    y : sequence
         input dimension of data points
    w : sequence or None, optional
        weights for spline fitting.  Must be positive.  If None (default),
        weights are all equal.
    bbox : sequence or None, optional
        2-sequence specifying the boundary of the approximation interval. If
        None (default), bbox=[x[0],x[-1]].
    k : int, optional
        Degree of the smoothing spline.  Must be <= 5.


    See Also
    --------
    UnivariateSpline : Superclass -- allows knots to be selected by a
        smoothing condition
    LSQUnivariateSpline : spline for which knots are user-selected
    splrep : An older, non object-oriented wrapping of FITPACK
    splev, sproot, splint, spalde
    BivariateSpline : A similar class for two-dimensional spline interpolation



    Examples
    --------
    >>> from numpy import linspace,exp
    >>> from numpy.random import randn
    >>> from scipy.interpolate import UnivariateSpline
    >>> x = linspace(-3,3,100)
    >>> y = exp(-x**2) + randn(100)/10
    >>> s = UnivariateSpline(x,y,s=1)
    >>> xs = linspace(-3,3,1000)
    >>> ys = s(xs)

    xs,ys is now a smoothed, super-sampled version of the noisy gaussian x,y

    """

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
    """
    One-dimensional spline with explicit internal knots.

    Fits a spline y=s(x) of degree `k` to the provided `x`,`y` data.  `t`
    specifies the internal knots of the spline

    Parameters
    ----------
    x : sequence
        input dimension of data points -- must be increasing
    y : sequence
        input dimension of data points
    t: sequence
        interior knots of the spline.  Must be in ascending order
        and bbox[0]<t[0]<...<t[-1]<bbox[-1]
    w : sequence or None, optional
        weights for spline fitting.  Must be positive.  If None (default),
        weights are all equal.
    bbox : sequence or None, optional
        2-sequence specifying the boundary of the approximation interval. If
        None (default), bbox=[x[0],x[-1]].
    k : int, optional
        Degree of the smoothing spline.  Must be <= 5.

    Raises
    ------

    ValueError
        If the interior knots do not satisfy the Schoenberg-Whitney conditions

    See Also
    --------
    UnivariateSpline : Superclass -- knots are specified by setting a
        smoothing condition
    InterpolatedUnivariateSpline : spline passing through all points
    splrep : An older, non object-oriented wrapping of FITPACK
    splev, sproot, splint, spalde
    BivariateSpline : A similar class for two-dimensional spline interpolation



    Examples
    --------
    >>> from numpy import linspace,exp
    >>> from numpy.random import randn
    >>> from scipy.interpolate import LSQUnivariateSpline
    >>> x = linspace(-3,3,100)
    >>> y = exp(-x**2) + randn(100)/10
    >>> t = [-1,0,1]
    >>> s = LSQUnivariateSpline(x,y,t)
    >>> xs = linspace(-3,3,1000)
    >>> ys = s(xs)

    xs,ys is now a smoothed, super-sampled version of the noisy gaussian x,y
    with knots [-3,-1,0,1,3]

    """

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
        if not alltrue(t[k+1:n-k]-t[k:n-k-1] > 0,axis=0):
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


class BivariateSpline(object):
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
    def __call__(self,x,y,mth='array'):
        """ Evaluate spline at positions x,y."""
        if mth=='array':
            tx,ty,c = self.tck[:3]
            kx,ky = self.degrees
            z,ier = dfitpack.bispev(tx,ty,c,kx,ky,x,y)
            assert ier==0,'Invalid input: ier='+`ier`
            return z
        raise NotImplementedError

    def ev(self, xi, yi):
        """
        Evaluate spline at points (x[i], y[i]), i=0,...,len(x)-1
        """
        tx,ty,c = self.tck[:3]
        kx,ky = self.degrees
        zi,ier = dfitpack.bispeu(tx,ty,c,kx,ky,xi,yi)
        assert ier==0, 'Invalid input: ier='+`ier`
        return zi

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
        return dfitpack.dblint(tx,ty,c,kx,ky,xa,xb,ya,yb)

class SmoothBivariateSpline(BivariateSpline):
    """ Smooth bivariate spline approximation.

    See also:

    bisplrep, bisplev - an older wrapping of FITPACK
    UnivariateSpline - a similar class for univariate spline interpolation
    LSQUnivariateSpline - to create a BivariateSpline using weighted
                          least-squares fitting
    """

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
        xb,xe,yb,ye = bbox
        nx,tx,ny,ty,c,fp,wrk1,ier = dfitpack.surfit_smth(x,y,z,w,
                                                         xb,xe,yb,ye,
                                                         kx,ky,s=s,
                                                         eps=eps,lwrk2=1)
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
    See also:

    bisplrep, bisplev - an older wrapping of FITPACK
    UnivariateSpline - a similar class for univariate spline interpolation
    SmoothUnivariateSpline - to create a BivariateSpline through the
                             given points
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
                       By default, bbox=[min(x,tx),max(x,tx),
                                         min(y,ty),max(y,ty)]
          kx,ky=3,3  - degrees of the bivariate spline.
          eps        - a threshold for determining the effective rank
                       of an over-determined linear system of
                       equations. 0 < eps < 1, default is 1e-16.
        """
        nx = 2*kx+2+len(tx)
        ny = 2*ky+2+len(ty)
        tx1 = zeros((nx,),float)
        ty1 = zeros((ny,),float)
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

class RectBivariateSpline(BivariateSpline):
    """ Bivariate spline approximation over a rectangular mesh.

    Can be used for both smoothing or interpolating data.

    See also:

    SmoothBivariateSpline - a smoothing bivariate spline for scattered data
    bisplrep, bisplev - an older wrapping of FITPACK
    UnivariateSpline - a similar class for univariate spline interpolation
    """

    def __init__(self, x, y, z,
                 bbox = [None]*4, kx=3, ky=3, s=0):
        """
        Input:
          x,y  - 1-d sequences of coordinates in strictly ascending order
            z  - 2-d array of data with shape (x.size,y.size)
        Optional input:
          bbox       - 4-sequence specifying the boundary of
                       the rectangular approximation domain.
                       By default, bbox=[min(x,tx),max(x,tx),
                                         min(y,ty),max(y,ty)]
          kx,ky=3,3  - degrees of the bivariate spline.
          s          - positive smoothing factor defined for
                       estimation condition:
                         sum((w[i]*(z[i]-s(x[i],y[i])))**2,axis=0) <= s
                       Default s=0 which is for interpolation
        """
        x,y = ravel(x),ravel(y)
        if not all(diff(x) > 0.0):
            raise TypeError,'x must be strictly increasing'
        if not all(diff(y) > 0.0):
            raise TypeError,'y must be strictly increasing'
        if not ((x.min() == x[0]) and (x.max() == x[-1])):
            raise TypeError, 'x must be strictly ascending'
        if not ((y.min() == y[0]) and (y.max() == y[-1])):
            raise TypeError, 'y must be strictly ascending'
        if not x.size == z.shape[0]:
            raise TypeError,\
                  'x dimension of z must have same number of elements as x'
        if not y.size == z.shape[1]:
            raise TypeError,\
                  'y dimension of z must have same number of elements as y'
        z = ravel(z)
        xb,xe,yb,ye = bbox
        nx,tx,ny,ty,c,fp,ier = dfitpack.regrid_smth(x,y,z,
                                                    xb,xe,yb,ye,
                                                    kx,ky,s)
        if ier in [0,-1,-2]: # normal return
            pass
        else:
            message = _surfit_messages.get(ier,'ier=%s' % (ier))
            warnings.warn(message)

        self.fp = fp
        self.tck = tx[:nx],ty[:ny],c[:(nx-kx-1)*(ny-ky-1)]
        self.degrees = kx,ky
