import numpy as N
import numpy.linalg as L
from model import LikelihoodModel, LikelihoodModelResults, ContrastResults

import utils
import scipy.linalg

class OLSModel(LikelihoodModel):
    
    """
    A simple ordinary least squares model.
    """

    def logL(self, b, Y, **extra):
        return -scipy.linalg.norm(self.whiten(Y) - N.dot(self.wdesign, b))**2 / 2.

    def __init__(self, design, **keywords):
        LikelihoodModel.__init__(self, **keywords)
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
        Estimate coefficients using lstsq, returning fitted values, Y
        and coefficients, but initialize is not called so no
        psuedo-inverse is calculated.
        """
            
        Z = self.whiten(Y)

        lfit = Results(L.lstsq(self.wdesign, Z)[0])
        lfit.predict = N.dot(self.design, lfit.beta)
        lfit.Y = Y

    def fit(self, Y, **keywords):
        """
        Full \'fit\' of the model including estimate of covariance matrix, (whitened)
        residuals and scale. 

        """
    
        Z = self.whiten(Y)
            

        lfit = Results(N.dot(self.calc_beta, Z),
                       normalized_cov_beta=self.normalized_cov_beta)

        lfit.df_resid = self.df_resid
        lfit.resid = Z - N.dot(self.design, lfit.beta)
        lfit.scale=N.add.reduce(lfit.resid**2) / lfit.df_resid
        lfit.predict = N.dot(self.design, lfit.beta)

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
        LikelihoodModel.__init__(self, **keywords)
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
        LikelihoodModel.__init__(self, **keywords)
        self.weights = weights
        self.initialize(design)

    def whiten(self, X):
        if X.ndim == 1:
            return X * N.sqrt(self.weights)
        elif X.ndim == 2:
            c = N.sqrt(self.weights)
            v = N.zeros(X.shape, N.float64)
            for i in range(X.shape[1]):
                v[:,i] = X[:,i] * c
            return v
    

class Results(LikelihoodModelResults):
    """
    This class summarizes the fit of a linear regression model.
    It handles the output of contrasts, estimates of covariance, etc.
    """

    def norm_resid(self):
        """
        Residuals, normalized to have unit length.

        Note: residuals are whitened residuals.

        """

        if not hasattr(self, 'resid'):
            raise ValueError, 'need normalized residuals to estimate standard deviation'

        sdd = utils.recipr(self.sd) / N.sqrt(self.df)
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

