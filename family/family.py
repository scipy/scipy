import numpy as N
from scipy.stats.models.family import links as L
from scipy.stats.models.family import varfuncs as V

class Family(object):

    """
    A class to model one-parameter exponential
    families.

    INPUTS:
       link      -- a Link instance
       variance  -- a variance function (models means as a function
                    of mean)

    """

    valid = [-N.inf, N.inf]

    tol = 1.0e-05
    links = []

    def _setlink(self, link):
        self._link = link
        if hasattr(self, "links"):
            if link not in self.links:
                raise ValueError, 'invalid link for family, should be in %s' % `self.links`

    def _getlink(self):
        return self._link

    link = property(_getlink, _setlink)

    def __init__(self, link, variance):

        self.link = link
        self.variance = variance

    def weights(self, mu):

        """
        Weights for IRLS step.

        w = 1 / (link'(mu)**2 * variance(mu))

        INPUTS:
           mu  -- mean parameter in exponential family

        OUTPUTS:
           w   -- weights used in WLS step of GLM/GAM fit

        """

        return 1. / (self.link.deriv(mu)**2 * self.variance(mu))

    def deviance(self, Y, mu, scale=1.):
        """
        Deviance of (Y,mu) pair. Deviance is usually defined
        as the difference

        DEV = (SUM_i -2 log Likelihood(Y_i,mu_i) + 2 log Likelihood(mu_i,mu_i)) / scale

        INPUTS:
           Y     -- response variable
           mu    -- mean parameter
           scale -- optional scale in denominator of deviance

        OUTPUTS: dev
           dev   -- DEV, as described aboce

        """

        return N.power(self.devresid(Y, mu), 2).sum() / scale

    def devresid(self, Y, mu):
        """
        The deviance residuals, defined as the residuals
        in the deviance.

        Without knowing the link, they default to Pearson residuals

        resid_P = (Y - mu) * sqrt(weight(mu))

        INPUTS:
           Y     -- response variable
           mu    -- mean parameter

        OUTPUTS: resid
           resid -- deviance residuals
        """

        return (Y - mu) * N.sqrt(self.weights(mu))

    def fitted(self, eta):
        """
        Fitted values based on linear predictors eta.

        INPUTS:
           eta  -- values of linear predictors, say,
                   X beta in a generalized linear model.

        OUTPUTS: mu
           mu   -- link.inverse(eta), mean parameter based on eta

        """
        return self.link.inverse(eta)

    def predict(self, mu):
        """
        Linear predictors based on given mu values.

        INPUTS:
           mu   -- mean parameter of one-parameter exponential family

        OUTPUTS: eta
           eta  -- link(mu), linear predictors, based on
                   mean parameters mu

        """
        return self.link(mu)

class Poisson(Family):

    """
    Poisson exponential family.

    INPUTS:
       link      -- a Link instance

    """

    links = [L.log, L.identity, L.sqrt]
    variance = V.mu
    valid = [0, N.inf]

    def __init__(self, link=L.log):
        self.variance = Poisson.variance
        self.link = link

    def devresid(self, Y, mu):
        """
        Poisson deviance residual

        INPUTS:
           Y     -- response variable
           mu    -- mean parameter

        OUTPUTS: resid
           resid -- deviance residuals

        """
        return N.sign(Y - mu) * N.sqrt(2 * Y * N.log(Y / mu) - 2 * (Y - mu))

class Gaussian(Family):

    """
    Gaussian exponential family.

    INPUTS:
       link      -- a Link instance

    """

    links = [L.log, L.identity, L.inverse]
    variance = V.constant

    def __init__(self, link=L.identity):
        self.variance = Gaussian.variance
        self.link = link

    def devresid(self, Y, mu, scale=1.):
        """
        Gaussian deviance residual

        INPUTS:
           Y     -- response variable
           mu    -- mean parameter
           scale -- optional scale in denominator (after taking sqrt)

        OUTPUTS: resid
           resid -- deviance residuals
        """

        return (Y - mu) / N.sqrt(self.variance(mu) * scale)

class Gamma(Family):

    """
    Gamma exponential family.

    INPUTS:
       link      -- a Link instance

    BUGS:
       no deviance residuals?

    """

    links = [L.log, L.identity, L.inverse]
    variance = V.mu_squared

    def __init__(self, link=L.identity):
        self.variance = Gamma.variance
        self.link = link


class Binomial(Family):

    """
    Binomial exponential family.

    INPUTS:
       link      -- a Link instance
       n         -- number of trials for Binomial
    """

    links = [L.logit, L.probit, L.cauchy, L.log, L.cloglog]
    variance = V.binary

    def __init__(self, link=L.logit, n=1):
        self.n = n
        self.variance = V.Binomial(n=self.n)
        self.link = link

    def devresid(self, Y, mu):
        """
        Binomial deviance residual

        INPUTS:
           Y     -- response variable
           mu    -- mean parameter

        OUTPUTS: resid
           resid -- deviance residuals

        """

        mu = self.link.clean(mu)
        return N.sign(Y - mu) * N.sqrt(-2 * (Y * N.log(mu / self.n) + (self.n - Y) * N.log(1 - mu / self.n)))

class InverseGaussian(Family):

    """
    InverseGaussian exponential family.

    INPUTS:
       link      -- a Link instance
       n         -- number of trials for Binomial

    """

    links = [L.inverse_squared, L.inverse, L.identity, L.log]
    variance = V.mu_cubed

    def __init__(self, link=L.identity):
        self.n = n
        self.variance = InverseGaussian.variance
        self.link = link
