import numpy as N
import numpy.linalg as L
from model import LikelihoodModel
import utils
import scipy.linalg

class OLSModel(LikelihoodModel):
    
    """
    A simple ordinary least squares model.
    """

    def logL(self, b, Y, **extra):
        return -scipy.linalg.norm(self.whiten(Y) - N.dot(self.wdesign, b))**2 / 2.

    def __init__(self, design, **keywords):
        Model.__init__(self, **keywords)
        self.initialize(design)

    def initialize(self, design, **keywords):
        self.design = design
        self.wdesign = self.whiten(design)
        self.calc_beta = L.pinv(self.wdesign)
        self.normalized_cov_beta = N.dot(self.calc_beta,
                                         N.transpose(self.calc_beta))
        self.df_resid = self.wdesign.shape[0] - utils.rank(self.design)

    def whiten(self, Y):
        return Y
    
    def est_coef(self, Y):
        """
        Estimate coefficients using lstsq, returning fitted values, Y and coefficients, but
        initialize is not called no psuedo-inverse is calculated.
        """
            
        lfit = Results()
        Z = self.whiten(Y)

        lfit.beta = L.lstsq(self.wdesign, Z)[0]
        lfit.predict = N.dot(self.design, lfit.beta)
        lfit.Y = Y

    def fit(self, Y, **keywords):
        """
        Full \'fit\' of the model including estimate of covariance matrix, (whitened)
        residuals and scale. 

        """
    
        Z = self.whiten(Y)
            
        lfit = Results()

        lfit.beta = N.dot(self.calc_beta, Z)
        lfit.normalized_cov_beta = self.normalized_cov_beta
        lfit.df_resid = self.df_resid
        lfit.predict = N.dot(self.design, lfit.beta)
        lfit.resid = Z - N.dot(self.design, lfit.beta)

        lfit.scale = N.add.reduce(lfit.resid**2) / lfit.df_resid

        lfit.Z = Z # just in case
        lfit.Y = Y
        
        return lfit

class ARModel(OLSModel):
    """
    A regression model with an AR(1) covariance structure.

    Eventually, this will be AR(p) -- all that is needed is to
    determine the self.whiten method from AR(p) parameters.
    """

    def __init__(self, design, rho=0, **keywords):
        Model.__init__(self, **keywords)
        self.rho = rho
        self.initialize(design)

    def whiten(self, X):
        factor = 1. / N.sqrt(1 - self.rho**2)
        return N.concatenate([[X[0]], (X[1:] - self.rho * X[0:-1]) * factor])

class WLSModel(ARModel):
    """

    A regression model with diagonal but non-identity covariance
    structure. The weights are proportional to the inverse of the
    variance of the observations.

    """

    def __init__(self, design, weights=1, **keywords):
        Model.__init__(self, **keywords)
        self.weights = weights
        self.initialize(design)

    def whiten(self, X):
        if X.ndim == 1:
            return X * N.sqrt(self.weights)
        elif X.ndim == 2:
            c = N.sqrt(self.weights)
            v = N.zeros(X.shape, N.Float)
            for i in range(X.shape[1]):
                v[:,i] = X[:,i] * c
            return v

class Results:
    """
    This class summarizes the fit of a linear regression model.
    It handles the output of contrasts, estimates of covariance, etc.
    """

    def t(self, column=None):
        """
        Return the t-statistic for a given column of the design matrix.

        Use Tcontrast for more complicated t-statistics.

        """
        if not hasattr(self, '_sd'):
            self.sd()
        if column is None:
            _t = N.zeros(_beta.shape, N.Float)
            for i in range(self.beta.shape[0]):
                _t[i] = _beta[i] * utils.recipr((self._sd * self.sqrt(self.normalized_cov_beta[i,i])))
        else:
            i = column
            _t = _beta[i] * utils.recipr((self._sd * self.sqrt(self.normalized_cov_beta[i,i])))
        return _t

    def sd(self):
        """
        Residual standard error: sqrt((resid**2).sum()) / df_resid.

        Note: residuals are whitened residuals.

        """
        if not hasattr(self, 'resid'):
            raise ValueError, 'need residuals to estimate standard deviation'
        self._sd = N.sqrt(self.scale)
        
    def norm_resid(self):
        """
        Residuals, normalized to have unit length.

        Note: residuals are whitened residuals.

        """

        if not hasattr(self, 'resid'):
            raise ValueError, 'need normalized residuals to estimate standard deviation'
        if not hasattr(self, '_sd'):
            self.sd()
        sdd = utils.recipr(self._sd) / N.sqrt(self.df)
        norm_resid = self.resid * N.multiply.outer(N.ones(Y.shape[0]), sdd)
        return norm_resid

    def predict(self, design):
        """
        Return fitted values from a design matrix.
        """

        return N.dot(design, self.beta)

    def Rsq(self, adjusted=False):
        """
        Return the R^2 value for each row of the response Y.
        """
        self.Ssq = N.std(self.Z)**2
        ratio = self.scale / self.Ssq
        if not adjusted: ratio *= ((Y.shape[0] - 1) / self.df_resid)
        return 1 - ratio

    def cov_beta(self, matrix=None, column=None, scale=None, other=None):
        """
        Returns the variance/covariance matrix of a linear contrast
        of the estimates of beta, multiplied by scale which
        will usually be an estimate of sigma^2.

        The covariance of
        interest is either specified as a (set of) column(s) or a matrix.
        """
        if scale is None:
            scale = 1.

        if column is not None:
            column = N.asarray(column)
            if column.shape == ():
                return self.normalized_cov_beta[column, column] * scale
            else:
                return self.normalized_cov_beta[column][:,column] * scale

        elif matrix is not None:
            if other is None:
                other = matrix
            tmp = N.dot(matrix, N.dot(self.normalized_cov_beta, N.transpose(other)))
            return tmp * scale

        if matrix is None and column is None:
            return self.normalized_cov_beta * scale

    def Tcontrast(self, matrix, t=True, sd=True, scale=None):
        """
        Compute a Tcontrast for a row vector matrix. To get the t-statistic
        for a single column, use the 't' method.
        """
        
        if not hasattr(self, '_sd'):
            self.sd()
        results = ContrastResults()

        results.effect = N.dot(matrix, self.beta)
        if sd:
            results.sd = N.sqrt(self.cov_beta(matrix=matrix)) * self._sd
        if t:
            results.t = results.effect * utils.recipr(results.sd)
        return results

    def Fcontrast(self, matrix, eff=True, t=True, sd=True, scale=None, invcov=None):
        """
        Compute an Fcontrast for a contrast matrix. 

        Here, matrix M is assumed to be non-singular. More precisely,

        M pX pX' M'

        is assumed invertible. Here, pX is the generalized inverse of the
        design matrix of the model. There can be problems in non-OLS models
        where the rank of the covariance of the noise is not full.

        See the contrast module to see how to specify contrasts.
        In particular, the matrices from these contrasts will always be
        non-singular in the sense above.

        """
        
        results = ContrastResults()
        cbeta = N.dot(matrix, self.beta)

        q = matrix.shape[0]
        if invcov is None:
            invcov = L.inv(self.cov_beta(matrix=matrix, scale=1.0))
        results.F = N.add.reduce(N.dot(invcov, cbeta) * cbeta, 0) * utils.recipr((q * self.scale))
        return results
    

def isestimable(C, D, pinv=None, warn=True):
    """
    From an q x p contrast matrix C and an n x p design matrix D, checks
    if the contrast C is estimable by looking at the rank of hstack([C,D]) and
    verifying it is the same as the rank of D.
    """

    if C.ndim == 1:
        C.shape = (C.shape[0], 1)
    new = N.hstack([C, D])
    if utils.rank(new) != utils.rank(D):
        return False
    return True
