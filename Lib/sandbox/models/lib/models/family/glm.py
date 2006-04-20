import numpy as N
import family
from models.regression import WLSModel

class GeneralizedLinearModel(WLSModel):

    niter = 10
    
    def __init__(self, design, family=family.Gaussian(), **keywords):
        self.family = family
        self.weights = 1
        self.initialize(design)

    def __iter__(self):
        self.iter = 0
        self.dev = N.inf
        return self

    def deviance(self, Y=None, results=None, scale = 1.):
        """
        Return (unnormalized) log-likelihood for glm.

        Note that self.scale is interpreted as a variance in OLSModel, so
        we divide the residuals by its sqrt.
        """
        if results is None:
            results = self.results
        if Y is None: Y = self.Y
        return self.family.deviance(Y, results.mu) / scale

    def next(self, results, Y):
        self.weights = self.family.weights(results.mu)
        self.initialize(self.design)
        Z = results.fitted + self.family.link.deriv(results.mu) * (Y - results.mu)
        newresults = WLSModel.fit(self, Z)
        newresults.mu = self.family.link.inverse(newresults.fitted)
        self.iter += 1
        return newresults

    def cont(self, results, tol=1.0e-05):
        """
        Continue iterating, or has convergence been obtained?
        """
        if self.iter >= GeneralizedLinearModel.niter:
            return False

        curdev = self.deviance(results=results)

        if N.fabs((self.dev - curdev) / curdev) < tol:
            return False
        self.dev = curdev
        
        return True

    def estimate_scale(self, Y=None, results=None):
        """
        Return Pearson\'s X^2 estimate of scale.
        """
        
        if results is None:
            results = self.results
        if Y is None: Y = self.Y
        resid = Y - results.mu
        return (N.power(resid, 2) / self.family.variance(results.mu)).sum() / results.df_resid
    
    def fit(self, Y, **keywords):

        self.Y = N.asarray(Y, N.Float)
        iter(self)
        self.results = WLSModel.fit(self, self.family.link(Y, initialize=True))
        self.results.mu = self.family.link.inverse(self.results.fitted)
        self.scale = self.results.scale = self.estimate_scale()
        
        while self.cont(self.results):
            self.results = self.next(self.results, Y)
            self.scale = self.results.scale = self.estimate_scale()

        return self.results
