import links as L
import varfuncs as V
import numpy as N

class Family:

    valid = [-N.inf, N.inf]

    def __init__(self, link, variance):
        self.link = link
        self.variance = variance

    def weights(self, mu):
        """
        Weights for IRLS step.
        """

        return 1. / (self.link.deriv(mu)**2 * self.variance(mu))

    def deviance(self, Y, mu, scale=1.):
        return N.power(self.devresid(Y, mu), 2).sum() / scale

    def devresid(self, Y, mu):
        return (Y - mu) * N.sqrt(self.weights(mu))

    def fitted(self, eta):
        """
        Fitted values based on linear predictors eta.
        """
        return self.link.inverse(eta)

    def predict(self, mu):
        """
        Linear predictors based on given mu values.
        """
        return self.link(mu)

class Poisson(Family):

    """
    Poisson exponential family in glm context.
    """

    links = [L.log, L.identity, L.sqrt]
    variance = V.mu
    valid = [0, N.inf]

    def __init__(self, link=L.log):
        if link not in Poisson.links:
            raise ValueError, 'invalid link for Poisson family'
        self.variance = Poisson.variance
        self.link = link

    def devresid(self, Y, mu):
        return N.sign(Y - mu) * N.sqrt(2 * Y * N.log(Y / mu) - 2 * (Y - mu))

class Gaussian(Family):

    """
    Gaussian exponential family in glm context.
    """

    links = [L.log, L.identity, L.inverse]
    variance = V.constant

    def __init__(self, link=L.identity):
        if link not in Gaussian.links:
            raise ValueError, 'invalid link for Gaussian family'
        self.variance = Gaussian.variance
        self.link = link

    def devresid(self, Y, mu, scale=1.):
        return (Y - mu) / N.sqrt(self.variance(mu) * scale)

class Gamma(Family):

    """
    Gaussian exponential family in glm context.
    """

    links = [L.log, L.identity, L.inverse]
    variance = V.mu_squared

    def __init__(self, link=L.identity):
        if link not in Gamma.links:
            raise ValueError, 'invalid link for Gamma family'
        self.variance = Gamma.variance
        self.link = link

    

class Binomial(Family):

    """
    Binomial exponential family in glm context.
    """

    links = [L.logit, L.probit, L.cauchy, L.log, L.cloglog]
    variance = V.binary

    def __init__(self, link=L.logit, n=1):
        if link not in Binomial.links:
            raise ValueError, 'invalid link for Binomial family'
        self.n = n
        self.variance = V.Binomial(n=self.n)
        self.link = link

    def devresid(self, Y, mu, scale=1.):
        mu = self.link.clean(mu)
        return N.sign(Y - mu) * N.sqrt(-2 * (Y * N.log(mu / self.n) + (self.n - Y) * N.log(1 - mu / self.n)))

class InverseGaussian(Family):

    """
    Gaussian exponential family in glm context.
    """

    links = [L.inverse_squared, L.inverse, L.identity, L.log]
    variance = V.mu_cubed

    def __init__(self, link=L.identity, n=1):
        if link not in InverseGaussian.links:
            raise ValueError, 'invalid link for InverseGaussian family'
        self.n = n
        self.variance = InverseGaussian.variance
        self.link = link

