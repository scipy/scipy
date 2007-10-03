"""
General linear models
--------------------

"""

import numpy as N
from scipy.stats.models import family
from scipy.stats.models.regression import WLSModel

class Model(WLSModel):

    niter = 10
    
    def __init__(self, design, family=family.Gaussian()):
        self.family = family
        super(Model, self).__init__(design, weights=1)

    def __iter__(self):                                                   
        self.iter = 0
        self.dev = N.inf
        return self

    def deviance(self, Y=None, results=None, scale = 1.):
        """
        Return (unnormalized) log-likelihood for GLM.

        Note that self.scale is interpreted as a variance in old_model, so
        we divide the residuals by its sqrt.
        """
        if results is None:
            results = self.results
        if Y is None:
            Y = self.Y
        return self.family.deviance(Y, results.mu) / scale

    def next(self):
        results = self.results
        Y = self.Y
        self.weights = self.family.weights(results.mu)
        self.initialize(self.design)
        Z = results.predict + self.family.link.deriv(results.mu) * (Y - results.mu)
        newresults = super(Model, self).fit(self, Z)
        newresults.Y = Y
        newresults.mu = self.family.link.inverse(newresults.predict)
        self.iter += 1
        return newresults

    def cont(self, tol=1.0e-05):
        """
        Continue iterating, or has convergence been obtained?
        """
        if self.iter >= Model.niter:
            return False

        curdev = self.deviance(results=self.results)

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
        if Y is None:
            Y = self.Y
        resid = Y - results.mu
        return ((N.power(resid, 2) / self.family.variance(results.mu)).sum()
                / results.df_resid)

    def fit(self, Y):
        self.Y = N.asarray(Y, N.float64)
        iter(self)
        self.results = super(Model, self).fit(
            self.family.link.initialize(Y))
        self.results.mu = self.family.link.inverse(self.results.predict)
        self.scale = self.results.scale = self.estimate_scale()
        
        while self.cont(self.results):
            self.results = self.next()
            self.scale = self.results.scale = self.estimate_scale()

        return self.results
