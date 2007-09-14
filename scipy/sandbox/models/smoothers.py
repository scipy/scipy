"""
This module contains scatterplot smoothers, that is classes
who generate a smooth fit of a set of (x,y) pairs.
"""

import numpy as N
import numpy.linalg as L

from scipy.linalg import solveh_banded
from scipy.optimize import golden

from scipy.sandbox.models import _bspline
from scipy.sandbox.models.bspline import bspline, _band2array


class poly_smoother:
    """
    Polynomial smoother up to a given order.
    Fit based on weighted least squares.

    The x values can be specified at instantiation or when called.
    """

    def df_fit(self):
        """
        Degrees of freedom used in the fit.
        """
        return self.order + 1

    def df_resid(self):
        """
        Residual degrees of freedom from last fit.
        """
        return self.N - self.order - 1

    def __init__(self, order, x=None):
        self.order = order
        self.coef = N.zeros((order+1,), N.float64)
        if x is not None:
            self.X = N.array([x**i for i in range(order+1)]).T

    def __call__(self, x=None):
        if x is not None:
            X = N.array([(x**i) for i in range(self.order+1)])
        else: X = self.X
        return N.squeeze(N.dot(X.T, self.coef)) 
    
    def fit(self, y, x=None, weights=None):
        self.N = y.shape[0]
        if weights is None:
            weights = 1
        _w = N.sqrt(weights)
        if x is None:
            if not hasattr(self, "X"):
                raise ValueError, "x needed to fit poly_smoother"
        else:
            self.X = N.array([(x**i) for i in range(self.order+1)])

        X = self.X * _w

        _y = y * _w
        self.coef = N.dot(L.pinv(X).T, _y)

class smoothing_spline(bspline):

    penmax = 30.

    def fit(self, y, x=None, weights=None, pen=0.):
        banded = True

        if x is None:
            x = self.tau[(self.M-1):-(self.M-1)] # internal knots

        if pen == 0.: # can't use cholesky for singular matrices
            banded = False
            
        if x.shape != y.shape:
            raise ValueError, 'x and y shape do not agree, by default x are the Bspline\'s internal knots'
        
        bt = self.basis(x)
        if pen >= self.penmax:
            pen = self.penmax

        if weights is None:
            weights = N.array(1.)

        wmean = weights.mean()
        _w = N.sqrt(weights / wmean)
        bt *= _w

        # throw out rows with zeros (this happens at boundary points!)

        mask = N.flatnonzero(1 - N.alltrue(N.equal(bt, 0), axis=0))

        bt = bt[:, mask]
        y = y[mask]

        self.df_total = y.shape[0]

        if bt.shape[1] != y.shape[0]:
            raise ValueError, "some x values are outside range of B-spline knots"
        bty = N.dot(bt, _w * y)
        self.N = y.shape[0]
        if not banded:
            self.btb = N.dot(bt, bt.T)
            _g = _band2array(self.g, lower=1, symmetric=True)
            self.coef, _, self.rank = L.lstsq(self.btb + pen*_g, bty)[0:3]
            self.rank = min(self.rank, self.btb.shape[0])
        else:
            self.btb = N.zeros(self.g.shape, N.float64)
            nband, nbasis = self.g.shape
            for i in range(nbasis):
                for k in range(min(nband, nbasis-i)):
                    self.btb[k, i] = (bt[i] * bt[i+k]).sum()

            bty.shape = (1, bty.shape[0])
            self.chol, self.coef = solveh_banded(self.btb + 
                                                 pen*self.g,
                                                 bty, lower=1)

        self.coef = N.squeeze(self.coef)
        self.resid = N.sqrt(wmean) * (y * _w - N.dot(self.coef, bt))
        self.pen = pen

    def gcv(self):
        """
        Generalized cross-validation score of current fit.
        """

        norm_resid = (self.resid**2).sum()
        return norm_resid / (self.df_total - self.trace())

    def df_resid(self):
        """
        self.N - self.trace()

        where self.N is the number of observations of last fit.
        """
        
        return self.N - self.trace()

    def df_fit(self):
        """
        = self.trace()

        How many degrees of freedom used in the fit? 
        """
        return self.trace()

    def trace(self):
        """
        Trace of the smoothing matrix S(pen)
        """

        if self.pen > 0:
            _invband = _bspline.invband(self.chol.copy())
            tr = _trace_symbanded(_invband, self.btb, lower=1)
            return tr
        else:
            return self.rank

class smoothing_spline_fixeddf(smoothing_spline):
    """
    Fit smoothing spline with approximately df degrees of freedom
    used in the fit, i.e. so that self.trace() is approximately df.

    In general, df must be greater than the dimension of the null space
    of the Gram inner product. For cubic smoothing splines, this means
    that df > 2.
    """

    target_df = 5

    def __init__(self, knots, order=4, coef=None, M=None, target_df=None):
        if target_df is not None:
            self.target_df = target_df
        bspline.__init__(self, knots, order=order, coef=coef, M=M)
        self.target_reached = False

    def fit(self, y, x=None, df=None, weights=None, tol=1.0e-03):

        df = df or self.target_df

        apen, bpen = 0, 1.0e-03
        olddf = y.shape[0] - self.m

        if not self.target_reached:
            while True:
                curpen = 0.5 * (apen + bpen)
                smoothing_spline.fit(self, y, x=x, weights=weights, pen=curpen)
                curdf = self.trace()
                if curdf > df:
                    apen, bpen = curpen, 2 * curpen
                else:
                    apen, bpen = apen, curpen
                    if apen >= self.penmax:
                        raise ValueError, "penalty too large, try setting penmax higher or decreasing df"
                if N.fabs(curdf - df) / df < tol: 
                    self.target_reached = True
                    break
        else:
            smoothing_spline.fit(self, y, x=x, weights=weights, pen=self.pen)

class smoothing_spline_gcv(smoothing_spline):

    """
    Fit smoothing spline trying to optimize GCV.

    Try to find a bracketing interval for scipy.optimize.golden
    based on bracket.

    It is probably best to use target_df instead, as it is 
    sometimes difficult to find a bracketing interval.

    """

    def fit(self, y, x=None, weights=None, tol=1.0e-03,
            bracket=(0,1.0e-03)):
    
        def _gcv(pen, y, x):
            smoothing_spline.fit(y, x=x, pen=N.exp(pen), weights=weights)
            a = self.gcv()
            return a

        a = golden(_gcv, args=(y,x), brack=(-100,20), tol=tol)

def _trace_symbanded(a,b, lower=0):
    """
    Compute the trace(a*b) for two upper or lower banded real symmetric matrices.
    """

    if lower:
        t = _zero_triband(a * b, lower=1)
        return t[0].sum() + 2 * t[1:].sum()
    else:
        t = _zero_triband(a * b, lower=0)
        return t[-1].sum() + 2 * t[:-1].sum()



def _zero_triband(a, lower=0):
    """
    Zero out unnecessary elements of a real symmetric banded matrix.
    """

    nrow, ncol = a.shape
    if lower: 
        for i in range(nrow): a[i,(ncol-i):] = 0.
    else:
        for i in range(nrow): a[i,0:i] = 0.
    return a
