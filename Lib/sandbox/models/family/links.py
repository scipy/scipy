import numpy as N
import scipy.stats

class Link:

    def initialize(self, Y):
        return N.asarray(Y).mean() * N.ones(Y.shape)

class Logit(Link):

    """
    The logit transform as a link function:

    g(x) = log(x / (1 - x))
    """

    tol = 1.0e-10

    def clean(self, p):
        return N.clip(p, Logit.tol, 1. - Logit.tol)

    def __call__(self, p):
        p = self.clean(p)
        return N.log(p / (1. - p))

    def inverse(self, z):
        t = N.exp(z)
        return t / (1. + t)

    def deriv(self, p):
        p = self.clean(p)
        return 1. / (p * (1 - p))

logit = Logit()

class Power(Link):

    """
    The power transform as a link function:

    g(x) = x**power

    """

    def __init__(self, power=1.):
        self.power = power

    def __call__(self, x):
        return N.power(x, self.power)

    def inverse(self, x):
        return N.power(x, 1. / self.power)

    def deriv(self, x):
        return self.power * N.power(x, self.power - 1)

inverse = Power(power=-1.)
inverse.__doc__ = """

The inverse transform as a link function:

g(x) = 1 / x
"""

sqrt = Power(power=0.5)
sqrt.__doc__ = """

The square-root transform as a link function:

g(x) = sqrt(x)
"""

inverse_squared = Power(power=-2.)
inverse_squared.__doc__ = """

The inverse squared transform as a link function:

g(x) = 1 / x**2
"""

identity = Power(power=1.)
identity.__doc__ = """

The identity transform as a link function:

g(x) = x
"""

class Log(Link):

    """
    The log transform as a link function:

    g(x) = log(x)

    """

    tol = 1.0e-10

    def clean(self, x):
        return N.clip(x, Logit.tol, N.inf)

    def __call__(self, x, **extra):
        x = self.clean(x)
        return N.log(x)

    def inverse(self, z):
        return N.exp(z)

    def deriv(self, x):
        x = self.clean(x)
        return 1. / x

log = Log()

class CDFLink(Logit):

    """
    The use the CDF of a scipy.stats distribution as a link function:

    g(x) = dbn.ppf(x)

    """

    def __init__(self, dbn=scipy.stats.norm):
        self.dbn = dbn

    def __call__(self, p):
        p = self.clean(p)
        return self.dbn.ppf(p)

    def inverse(self, z):
        return self.dbn.cdf(z)

    def deriv(self, p):
        p = self.clean(p)
        return 1. / self.dbn.pdf(self(p))

probit = CDFLink()
probit.__doc__ = """

The probit (standard normal CDF) transform as a link function:

g(x) = scipy.stats.norm.ppf(x)

"""

cauchy = CDFLink(dbn=scipy.stats.cauchy)
cauchy.__doc__ = """

The Cauchy (standard Cauchy CDF) transform as a link function:

g(x) = scipy.stats.cauchy.ppf(x)

"""

class CLogLog(Logit):

    """
    The complementary log-log transform as a link function:

    g(x) = log(-log(x))

    """

    def __call__(self, p):
        p = self.clean(p)
        return N.log(-N.log(p))

    def inverse(self, z):
        return N.exp(-N.exp(z))

    def deriv(self, p):
        p = self.clean(p)
        return -1. / (N.log(p) * p)

cloglog = CLogLog()

