"""
This module implements some standard regression models: OLS and WLS
models, as well as an AR(p) regression model.

Models are specified with a design matrix and are fit using their 
'fit' method. 

Subclasses that have more complicated covariance matrices
should write over the 'whiten' method as the fit method
prewhitens the response by calling 'whiten'.

General reference for regression models:

'Introduction to Linear Regression Analysis', Douglas C. Montgomery,
    Elizabeth A. Peck, G. Geoffrey Vining. Wiley, 2006.

"""

__docformat__ = 'restructuredtext en'

import numpy as N
import numpy.linalg as L
from scipy.linalg import norm, toeplitz

from scipy.stats.models.model import likelihood_model, \
     likelihood_model_results
from scipy.stats.models import utils

class ols_model(likelihood_model):
    """
    A simple ordinary least squares model.

    Examples
    --------
    >>> import numpy as N
    >>> 
    >>> from scipy.stats.models.formula import term, I
    >>> from scipy.stats.models.regression import ols_model
    >>> 
    >>> data={'Y':[1,3,4,5,2,3,4],
    ...       'X':range(1,8)}
    >>> f = term("X") + I
    >>> f.namespace = data
    >>> 
    >>> model = ols_model(f.design())
    >>> results = model.fit(data['Y'])
    >>> 
    >>> results.beta
    array([ 0.25      ,  2.14285714])
    >>> results.t()
    array([ 0.98019606,  1.87867287])
    >>> print results.Tcontrast([0,1])
    <T contrast: effect=2.14285714286, sd=1.14062281591, t=1.87867287326, df_denom=5>
    >>> print results.Fcontrast(N.identity(2))
    <F contrast: F=19.4607843137, df_denom=5, df_num=2>
    """

    def logL(self, b, Y):
        return -norm(self.whiten(Y) - N.dot(self.wdesign, b))**2 / 2.

    def __init__(self, design):
        """
        Create a `ols_model` from a design.

        :Parameters:
            design : TODO
                TODO
        """
        super(ols_model, self).__init__()
        self.initialize(design)

    def initialize(self, design):
        """
        Set design for model, prewhitening design matrix and precomputing
        covariance of coefficients (up to scale factor in front).

        :Parameters:
            design : TODO
                TODO
        """

        self.design = design
        self.wdesign = self.whiten(design)
        self.calc_beta = L.pinv(self.wdesign)
        self.normalized_cov_beta = N.dot(self.calc_beta,
                                         N.transpose(self.calc_beta))
        self.df_resid = self.wdesign.shape[0] - utils.rank(self.design)

    def whiten(self, Y):
        """
        OLS model whitener does nothing: returns Y.
        """
        return Y
    
    def est_coef(self, Y):
        """
        Estimate coefficients using lstsq, returning fitted values, Y
        and coefficients, but initialize is not called so no
        psuedo-inverse is calculated.
        """
        Z = self.whiten(Y)

        lfit = regression_results(L.lstsq(self.wdesign, Z)[0], Y)
        lfit.predict = N.dot(self.design, lfit.beta)


    def fit(self, Y):
        """
        Full fit of the model including estimate of covariance matrix, 
        (whitened) residuals and scale. 

        """
        Z = self.whiten(Y)

        lfit = regression_results(N.dot(self.calc_beta, Z), Y,
                       normalized_cov_beta=self.normalized_cov_beta)

        lfit.df_resid = self.df_resid
        lfit.predict = N.dot(self.design, lfit.beta)
        lfit.resid = Z - N.dot(self.wdesign, lfit.beta)
        lfit.scale = N.add.reduce(lfit.resid**2) / lfit.df_resid

        lfit.Z = Z 
        
        return lfit

class ar_model(ols_model):
    """
    A regression model with an AR(p) covariance structure.

    The linear autoregressive process of order p--AR(p)--is defined as:
        TODO

    Examples
    --------
    >>> import numpy as N
    >>> import numpy.random as R
    >>> 
    >>> from scipy.stats.models.formula import term, I
    >>> from scipy.stats.models.regression import ar_model
    >>> 
    >>> data={'Y':[1,3,4,5,8,10,9],
    ...       'X':range(1,8)}
    >>> f = term("X") + I
    >>> f.namespace = data
    >>> 
    >>> model = ar_model(f.design(), 2)
    >>> for i in range(6):
    ...     results = model.fit(data['Y'])
    ...     print "AR coefficients:", model.rho
    ...     rho, sigma = model.yule_walker(data["Y"] - results.predict)
    ...     model = ar_model(model.design, rho)
    ... 
    AR coefficients: [ 0.  0.]
    AR coefficients: [-0.52571491 -0.84496178]
    AR coefficients: [-0.620642   -0.88654567]
    AR coefficients: [-0.61887622 -0.88137957]
    AR coefficients: [-0.61894058 -0.88152761]
    AR coefficients: [-0.61893842 -0.88152263]
    >>> results.beta
    array([ 1.58747943, -0.56145497])
    >>> results.t()
    array([ 30.796394  ,  -2.66543144])
    >>> print results.Tcontrast([0,1])
    <T contrast: effect=-0.561454972239, sd=0.210643186553, t=-2.66543144085, df_denom=5>
    >>> print results.Fcontrast(N.identity(2))
    <F contrast: F=2762.42812716, df_denom=5, df_num=2>
    >>>
    >>> model.rho = N.array([0,0])
    >>> model.iterative_fit(data['Y'], niter=3)
    >>> print model.rho
    [-0.61887622 -0.88137957]
    """
    def __init__(self, design, rho):
        if type(rho) is type(1):
            self.order = rho
            self.rho = N.zeros(self.order, N.float64)
        else:
            self.rho = N.squeeze(N.asarray(rho))
            if len(self.rho.shape) not in [0,1]:
                raise ValueError, "AR parameters must be a scalar or a vector"
            if self.rho.shape == ():
                self.rho.shape = (1,)
            self.order = self.rho.shape[0]
        super(ar_model, self).__init__(design)

    def iterative_fit(self, Y, niter=3):
        """
        Perform an iterative two-stage procedure to estimate AR(p)
        parameters and regression coefficients simultaneously.

        :Parameters:
            Y : TODO
                TODO
            niter : ``integer``
                the number of iterations
        """
        for i in range(niter):
            self.initialize(self.design)
            results = self.fit(Y)
            self.rho, _ = self.yule_walker(Y - results.predict)

    def whiten(self, X):
        """
        Whiten a series of columns according to an AR(p)
        covariance structure.

        :Parameters:
            X : TODO
                TODO
        """
        X = N.asarray(X, N.float64)
        _X = X.copy()
        for i in range(self.order):
            _X[(i+1):] = _X[(i+1):] - self.rho[i] * X[0:-(i+1)]
        return _X

    def yule_walker(self, X, method="unbiased", df=None):
        """
        Estimate AR(p) parameters from a sequence X using Yule-Walker equation.

        unbiased or maximum-likelihood estimator (mle)

        See, for example:

        http://en.wikipedia.org/wiki/Autoregressive_moving_average_model

        :Parameters:
            X : TODO
                TODO
            method : ``string``
                Method can be "unbiased" or "mle" and this determines
                denominator in estimate of autocorrelation function (ACF)
                at lag k. If "mle", the denominator is n=r.shape[0], if
                "unbiased" the denominator is n-k.
            df : ``integer``
                Specifies the degrees of freedom. If df is supplied,
                then it is assumed the X has df degrees of
                freedom rather than n.
        """
        
        method = str(method).lower()
        if method not in ["unbiased", "mle"]:
            raise ValueError, "ACF estimation method must be 'unbiased' \
            or 'MLE'"
        X = N.asarray(X, N.float64)
        X -= X.mean()
        n = df or X.shape[0]

        if method == "unbiased":
            denom = lambda k: n - k
        else:
            denom = lambda k: n

        if len(X.shape) != 1:
            raise ValueError, "expecting a vector to estimate AR parameters"
        r = N.zeros(self.order+1, N.float64)
        r[0] = (X**2).sum() / denom(0)
        for k in range(1,self.order+1):
            r[k] = (X[0:-k]*X[k:]).sum() / denom(k)
        R = toeplitz(r[:-1])

        rho = L.solve(R, r[1:])
        sigmasq = r[0] - (r[1:]*rho).sum()
        return rho, N.sqrt(sigmasq)

class wls_model(ols_model):
    """
    A regression model with diagonal but non-identity covariance
    structure. The weights are presumed to be
    (proportional to the) inverse of the
    variance of the observations.

    >>> import numpy as N
    >>> 
    >>> from scipy.stats.models.formula import term, I
    >>> from scipy.stats.models.regression import wls_model
    >>> 
    >>> data={'Y':[1,3,4,5,2,3,4],
    ...       'X':range(1,8)}
    >>> f = term("X") + I
    >>> f.namespace = data
    >>> 
    >>> model = wls_model(f.design(), weights=range(1,8))
    >>> results = model.fit(data['Y'])
    >>> 
    >>> results.beta
    array([ 0.0952381 ,  2.91666667])
    >>> results.t()
    array([ 0.35684428,  2.0652652 ])
    >>> print results.Tcontrast([0,1])
    <T contrast: effect=2.91666666667, sd=1.41224801095, t=2.06526519708, df_denom=5>
    >>> print results.Fcontrast(N.identity(2))
    <F contrast: F=26.9986072423, df_denom=5, df_num=2>
    """
    def __init__(self, design, weights=1):
        weights = N.array(weights)
        if weights.shape == (): # scalar
            self.weights = weights
        else: 
            design_rows = design.shape[0]
            if not(weights.shape[0] == design_rows and
                   weights.size == design_rows) :
                raise ValueError(
                    'Weights must be scalar or same length as design')
            self.weights = weights.reshape(design_rows)
        super(wls_model, self).__init__(design)

    def whiten(self, X):
        """
        Whitener for WLS model, multiplies by sqrt(self.weights)
        """
        X = N.asarray(X, N.float64)

        if X.ndim == 1:
            return X * N.sqrt(self.weights)
        elif X.ndim == 2:
            c = N.sqrt(self.weights)
            v = N.zeros(X.shape, N.float64)
            for i in range(X.shape[1]):
                v[:,i] = X[:,i] * c
            return v
    
class regression_results(likelihood_model_results):
    """
    This class summarizes the fit of a linear regression model.

    It handles the output of contrasts, estimates of covariance, etc.
    """

    def __init__(self, beta, Y, normalized_cov_beta=None, scale=1.):
        super(regression_results, self).__init__(beta,
                                                 normalized_cov_beta,
                                                 scale)
        self.Y = Y

    def norm_resid(self):
        """
        Residuals, normalized to have unit length.

        Note: residuals are whitened residuals.
        """
        if not hasattr(self, 'resid'):
            raise ValueError, 'need normalized residuals to estimate standard deviation'

        sdd = utils.recipr(self.sd) / N.sqrt(self.df)
        return  self.resid * N.multiply.outer(N.ones(self.Y.shape[0]), sdd)


    def predictors(self, design):
        """
        Return linear predictor values from a design matrix.
        """
        return N.dot(design, self.beta)

    def Rsq(self, adjusted=False):
        """
        Return the R^2 value for each row of the response Y.
        """
        self.Ssq = N.std(self.Z,axis=0)**2
        ratio = self.scale / self.Ssq
        if not adjusted: ratio *= ((self.Y.shape[0] - 1) / self.df_resid)
        return 1 - ratio

def isestimable(C, D):
    """
    From an q x p contrast matrix C and an n x p design matrix D, checks
    if the contrast C is estimable by looking at the rank of vstack([C,D]) and
    verifying it is the same as the rank of D.
    """
    if C.ndim == 1:
        C.shape = (C.shape[0], 1)
    new = N.vstack([C, D])
    if utils.rank(new) != utils.rank(D):
        return False
    return True

