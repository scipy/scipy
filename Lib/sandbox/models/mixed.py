import numpy as N
import numpy.linalg as L
from scipy.sandbox.models.formula import Formula, I

class Unit:

    """
    Individual experimental unit for 
    EM implementation of (repeated measures)
    mixed effects model.

    \'Maximum Likelihood Computations with Repeated Measures:
    Application of the EM Algorithm\'

    Nan Laird; Nicholas Lange; Daniel Stram 

    Journal of the American Statistical Association,
    Vol. 82, No. 397. (Mar., 1987), pp. 97-105. 
    """

    def __getitem__(self, item):
        return self.dict[item]

    def __setitem__(self, item, value):
        self.dict[item] = value

    def __init__(self, dict):
        self.dict = dict

    def __call__(self, formula, **extra):
        """
        Return the corresponding design matrix from formula,
        perform a check whether formula just has an intercept in it, in
        which case the number of rows must be computed.
        """
        if hasattr(self, 'n') and not extra.has_key('nrow'):
            extra['nrow'] = self.n
        return formula(namespace=self.dict, **extra)

    def design(self, formula, **extra):
        v = N.transpose(self(formula, **extra))
        self.n = v.shape[0]
        return v

    def compute_S(self, D, sigma):
        """
        Display (3.3) from Laird, Lange, Stram (see help(Unit))
        """ 

        self.S = (N.identity(self.n) * sigma**2 +
                  N.dot(self.Z, N.dot(D, N.transpose(self.Z))))

    def compute_W(self):
        """
        Display (3.2) from Laird, Lange, Stram (see help(Unit))
        """ 
        self.W = L.inv(self.S)

    def compute_P(self, Sinv):
        """
        Display (3.10) from Laird, Lange, Stram (see help(Unit))
        """ 
        t = N.dot(self.W, self.X)
        self.P = self.W - N.dot(N.dot(t, Sinv), N.transpose(t))

    def compute_r(self, alpha):
        """
        Display (3.5) from Laird, Lange, Stram (see help(Unit))
        """ 
        self.r = self.Y - N.dot(self.X, alpha)

    def compute_b(self, D):
        """
        Display (3.4) from Laird, Lange, Stram (see help(Unit))
        """ 
        self.b = N.dot(D, N.dot(N.dot(N.transpose(self.Z), self.W), self.r))

    def fit(self, a, D, sigma):
        """
        Compute unit specific parameters in
        Laird, Lange, Stram (see help(Unit)).

        Displays (3.2)-(3.5).
        """

        self.compute_S(D, sigma)
        self.compute_W()
        self.compute_r(a)
        self.compute_b(D)

    def compute_xtwy(self):
        """
        Utility function to compute X^tWY for Unit instance.
        """
        return N.dot(N.dot(self.W, self.Y), self.X)

    def compute_xtwx(self):
        """
        Utility function to compute X^tWX for Unit instance.
        """
        return N.dot(N.dot(N.transpose(self.X), self.W), self.X)

    def cov_random(self, D, Sinv=None):
        """
        Approximate covariance of estimates of random effects. Just after
        Display (3.10) in Laird, Lange, Stram (see help(Unit)).
        """
        if Sinv is not None:
            self.compute_P(Sinv)
        t = N.dot(self.Z, D)
        return D - N.dot(N.dot(N.transpose(t), self.P), t)

    def logL(self, a, ML=False):
        """
        Individual contributions to the log-likelihood, tries to return REML
        contribution by default though this requires estimated
        fixed effect a to be passed as an argument.
        """
        
        if ML:
            return (N.log(L.det(self.W)) - (self.r * N.dot(self.W, self.r)).sum()) / 2.
        else:
            if a is None:
                raise ValueError, 'need fixed effect a for REML contribution to log-likelihood'
            r = self.Y - N.dot(self.X, a)
            return (N.log(L.det(self.W)) - (r * N.dot(self.W, r)).sum()) / 2.

    def deviance(self, ML=False):
        return - 2 * self.logL(ML=ML)

class Mixed:

    """
    Model for 
    EM implementation of (repeated measures)
    mixed effects model.

    \'Maximum Likelihood Computations with Repeated Measures:
    Application of the EM Algorithm\'

    Nan Laird; Nicholas Lange; Daniel Stram 

    Journal of the American Statistical Association,
    Vol. 82, No. 397. (Mar., 1987), pp. 97-105. 
    """

    def __init__(self, units, response, fixed=I, random=I):
        self.units = units
        self.m = len(self.units)
        
        self.fixed = Formula(fixed)
        self.random = Formula(random)
        self.response = Formula(response)

        self.N = 0
        for unit in self.units:
            unit.Y = N.squeeze(unit.design(self.response))
            unit.X = unit.design(self.fixed)
            unit.Z = unit.design(self.random)
            self.N += unit.X.shape[0]

        # Determine size of fixed effects

        d = self.units[0].design(self.fixed)
        self.p = d.shape[1]  # d.shape = p 
        self.a = N.zeros(self.p, N.float64)

        # Determine size of D, and sensible initial estimates
        # of sigma and D
        
        d = self.units[0].design(self.random)
        self.q = d.shape[1]  # d.shape = q 

        self.D = N.zeros((self.q,)*2, N.float64)
        self.sigma = 1.

        self.dev = N.inf

    def compute_a(self):
        """
        Display (3.1) of
        Laird, Lange, Stram (see help(Mixed)).

        """

        S = 0
        Y = 0
        for unit in self.units:
            unit.fit(self.a, self.D, self.sigma)
            S += unit.compute_xtwx()
            Y += unit.compute_xtwy()
        self.Sinv = L.pinv(S)
        self.a = N.dot(self.Sinv, Y)
            
    def compute_sigma(self, ML=False):
        """
        Estimate sigma. If ML is True, return the ML estimate of sigma,
        else return the REML estimate.

        If ML, this is (3.6) in Laird, Lange, Stram (see help(Mixed)),
        otherwise it corresponds to (3.8).

        """
        sigmasq = 0.
        for unit in self.units:
            if ML:
                W = unit.W
            else:
                unit.compute_P(self.Sinv)
                W = unit.P
            t = unit.r - N.dot(unit.Z, unit.b)
            sigmasq += N.power(t, 2).sum()
            sigmasq += self.sigma**2 * N.trace(N.identity(unit.n) -
                                               self.sigma**2 * W)
        self.sigma = N.sqrt(sigmasq / self.N)

    def compute_D(self, ML=False):
        """
        Estimate random effects covariance D.
        If ML is True, return the ML estimate of sigma,
        else return the REML estimate.

        If ML, this is (3.7) in Laird, Lange, Stram (see help(Mixed)),
        otherwise it corresponds to (3.9).

        """
        D = 0.
        for unit in self.units:
            if ML:
                W = unit.W
            else:
                unit.compute_P(self.Sinv)
                W = unit.P
            D += N.multiply.outer(unit.b, unit.b)
            t = N.dot(unit.Z, self.D)
            D += self.D - N.dot(N.dot(N.transpose(t), W), t)

        self.D = D / self.m

    def cov_fixed(self):
        """
        Approximate covariance of estimates of fixed effects. Just after
        Display (3.10) in Laird, Lange, Stram (see help(Mixed)).
        """
        return self.Sinv

    def deviance(self, ML=False):
        return - 2 * self.logL(ML=ML)

    def logL(self, ML=False):
        """
        Return log-likelihood, REML by default.
        """
        logL = 0.

        for unit in self.units:
            logL += unit.logL(a=self.a, ML=ML)
        if not ML:
            logL += N.log(L.det(self.Sinv)) / 2
        return logL

    def initialize(self):
        S = 0
        Y = 0
        for unit in self.units:
            S += N.dot(N.transpose(unit.X), unit.X)
            Y += N.dot(N.transpose(unit.X), unit.Y)

        self.a = L.lstsq(S, Y)[0]

        D = 0
        t = 0
        sigmasq = 0
        for unit in self.units:
            unit.r = unit.Y - N.dot(unit.X, self.a)
            if self.q > 1:
                unit.b = L.lstsq(unit.Z, unit.r)[0]
            else:
                Z = unit.Z.reshape((unit.Z.shape[0], 1))
                unit.b = L.lstsq(Z, unit.r)[0]

            sigmasq += (N.power(unit.Y, 2).sum() -
                        (self.a * N.dot(N.transpose(unit.X), unit.Y)).sum() -
                        (unit.b * N.dot(N.transpose(unit.Z), unit.r)).sum())
            D += N.multiply.outer(unit.b, unit.b)
            t += L.pinv(N.dot(N.transpose(unit.Z), unit.Z))

        sigmasq /= (self.N - (self.m - 1) * self.q - self.p)
        self.sigma = N.sqrt(sigmasq)
        self.D = (D - sigmasq * t) / self.m

    def cont(self, ML=False, tol=1.0e-05):

        self.dev, old = self.deviance(ML=ML), self.dev
        if N.fabs((self.dev - old)) * self.dev < tol:
            return False
        return True

    def fit(self, niter=100, ML=False):

        for i in range(niter):
            self.compute_a()
            self.compute_sigma(ML=ML)
            self.compute_D(ML=ML)
            if not self.cont(ML=ML):
                break
            

if __name__ == '__main__':
    import numpy.random as R

    nsubj = 400
    units  = []
    
    n = 3

    import formula
    fixed = formula.Term('f')
    random = formula.Term('r')
    response = formula.Term('y')

    for i in range(nsubj):
        d = R.standard_normal()
        X = R.standard_normal((10,n))
        Z = X[0:2]
        Y = R.standard_normal((n,)) + d * 4
        units.append(Unit({'f':X, 'r':Z, 'y':Y}))

    m = Mixed(units, response)#, fixed, random)
    m.initialize()
    m.fit()

## a = Unit()
## a['x'] = N.array([2,3])
## a['y'] = N.array([3,4])

## x = Term('x')
## y = Term('y')

## fixed = x + y + x * y
## random = Formula(x)

## a.X = a.design(fixed)
## a.Z = a.design(random)

## print help(a.compute_S)


