"""
fitpack --- curve and surface fitting with splines

fitpack is based on a collection of Fortran routines DIERCKX
by P. Dierckx (see http://www.netlib.org/dierckx/) transformed
to double routines by Pearu Peterson.
"""
# Created by Pearu Peterson, June 2003

__all__ = ['LSQBivariateSpline','SmoothBivariateSpline']

import warnings
from scipy_base import zeros, Float
import dfitpack

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
        approximation.
        """
        return self.fp
    def get_knots(self):
        """ Return knots positions of the spline."""
        return self.tck[:2]
    def get_coeffs(self):
        """ Return spline coefficients."""
        return self.tck[2]
    def __call__(self,x,y,mth='array'):
        """ Evaluate spline at positions x,y."""

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
