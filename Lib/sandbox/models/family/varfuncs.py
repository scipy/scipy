import numpy as N

class VarianceFunction:
    """
    Variance function that relates the variance of a random variable
    to its mean. Defaults to 1.

    """

    def __call__(self, mu):
        return N.ones(mu.shape, N.float64)

constant = VarianceFunction()

class Power:

    """
    Variance function:

    V(mu) = fabs(mu)**power
    """

    def __init__(self, power=1.):
        self.power = power

    def __call__(self, mu):
        return N.power(N.fabs(mu), self.power)

class Binomial:

    """
    Binomial variance function

    p = mu / n; V(mu) = p * (1 - p) * n
    """

    tol = 1.0e-10

    def __init__(self, n=1):
        self.n = n

    def clean(self, p):
        return N.clip(p, Binomial.tol, 1 - Binomial.tol)

    def __call__(self, mu):
        p = self.clean(mu / self.n)
        return p * (1 - p) * self.n

mu = Power()
mu_squared = Power(power=2)
mu_cubed = Power(power=3)
binary = Binomial()
    
