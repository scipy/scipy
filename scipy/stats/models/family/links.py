import numpy as N
import scipy.stats

class Link:

    """
    A generic link function for one-parameter exponential
    family, with call, inverse and deriv methods.

    """

    def initialize(self, Y):
        return N.asarray(Y).mean() * N.ones(Y.shape)

    def __call__(self, p):
        return NotImplementedError

    def inverse(self, z):
        return NotImplementedError

    def deriv(self, p):
        return NotImplementedError


class Logit(Link):

    """
    The logit transform as a link function:

    g'(x) = 1 / (x * (1 - x))
    g^(-1)(x) = exp(x)/(1 + exp(x))

    """

    tol = 1.0e-10

    def clean(self, p):
        """
        Clip logistic values to range (tol, 1-tol)

        INPUTS:
           p     -- probabilities

        OUTPUTS: pclip
           pclip -- clipped probabilities
        """

        return N.clip(p, Logit.tol, 1. - Logit.tol)

    def __call__(self, p):
        """
        Logit transform

        g(p) = log(p / (1 - p))

        INPUTS:
           p   -- probabilities

        OUTPUTS: z
           z   -- logit transform of p

        """

        p = self.clean(p)
        return N.log(p / (1. - p))

    def inverse(self, z):
        """
        Inverse logit transform

        h(z) = exp(z)/(1+exp(z))

        INPUTS:
           z -- logit transform of p

        OUTPUTS: p
           p   -- probabilities

        """
        t = N.exp(z)
        return t / (1. + t)

    def deriv(self, p):

        """
        Derivative of logit transform

        g(p) = 1 / (p * (1 - p))

        INPUTS:
           p   -- probabilities

        OUTPUTS: y
           y   -- derivative of logit transform of p

        """
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
        """
        Power transform

        g(x) = x**self.power

        INPUTS:
           x   -- mean parameters

        OUTPUTS: z
           z   -- power transform of x

        """

        return N.power(x, self.power)

    def inverse(self, z):
        """
        Inverse of power transform

        g(x) = x**(1/self.power)

        INPUTS:
           z   -- linear predictors in GLM

        OUTPUTS: x
           x   -- mean parameters

        """
        return N.power(z, 1. / self.power)

    def deriv(self, x):
        """
        Derivative of power transform

        g(x) = self.power * x**(self.power - 1)

        INPUTS:
           x   -- mean parameters

        OUTPUTS: z
           z   -- derivative of power transform of x

        """

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
        """
        Log transform

        g(x) = log(x)

        INPUTS:
           x   -- mean parameters

        OUTPUTS: z
           z   -- log(x)

        """
        x = self.clean(x)
        return N.log(x)

    def inverse(self, z):
        """
        Inverse of log transform

        g(x) = exp(x)

        INPUTS:
           z   -- linear predictors in GLM

        OUTPUTS: x
           x   -- exp(z)

        """
        return N.exp(z)

    def deriv(self, x):
        """
        Derivative of log transform

        g(x) = 1/x

        INPUTS:
           x   -- mean parameters

        OUTPUTS: z
           z   -- derivative of log transform of x

        """

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
        """
        CDF link

        g(p) = self.dbn.pdf(p)

        INPUTS:
           p   -- mean parameters

        OUTPUTS: z
           z   -- derivative of CDF transform of p

        """
        p = self.clean(p)
        return self.dbn.ppf(p)

    def inverse(self, z):
        """
        Derivative of CDF link

        g(z) = self.dbn.cdf(z)

        INPUTS:
           z   -- linear predictors in GLM

        OUTPUTS: p
           p   -- inverse of CDF link of z

        """
        return self.dbn.cdf(z)

    def deriv(self, p):
        """
        Derivative of CDF link

        g(p) = 1/self.dbn.pdf(self.dbn.ppf(p))

        INPUTS:
           x   -- mean parameters

        OUTPUTS: z
           z   -- derivative of CDF transform of x

        """
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
        """
        C-Log-Log transform

        g(p) = log(-log(p))

        INPUTS:
           p   -- mean parameters

        OUTPUTS: z
           z   -- log(-log(p))

        """
        p = self.clean(p)
        return N.log(-N.log(p))

    def inverse(self, z):
        """
        Inverse of C-Log-Log transform

        g(z) = exp(-exp(z))

        INPUTS:
           z   -- linear predictor scale

        OUTPUTS: p
           p   -- mean parameters

        """
        return N.exp(-N.exp(z))

    def deriv(self, p):
        """
        Derivatve of C-Log-Log transform

        g(p) = - 1 / (log(p) * p)

        INPUTS:
           p   -- mean parameters

        OUTPUTS: z
           z   --  - 1 / (log(p) * p)

        """
        p = self.clean(p)
        return -1. / (N.log(p) * p)

cloglog = CLogLog()
