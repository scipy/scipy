import scipy.optimize
import numpy as N
from numpy.linalg import inv

from utils import recipr

class Model:

    """
    A (predictive) statistical model. The class Model itself does nothing
    but lays out the methods expected of any subclass.
    """

    def __init__(self, **keywords):
        pass

    def initialize(self, **keywords):
        """
        Initialize (possibly re-initialize) a Model instance. For
        instance, the design matrix of a linear model may change
        and some things must be recomputed.
        """
        raise NotImplementedError

    def fit(self, **keywords): 
        """
        Fit a model to data.
        """
        raise NotImplementedError

    def predict(self, **keywords):
        """
        After a model has been fit, results are (assumed to be) stored
        in self.results, which itself should have a predict method.
        """
        self.results.predict(**keywords) 

    def view(self, **keywords):
        """
        View results of a model.
        """
        raise NotImplementedError

class LikelihoodModel(Model):

    def logL(self, theta, **extra):
        """
        Log-likelihood of model.
        """
        raise NotImplementedError

    def score(self, theta, **extra):
        """
        Score function of model = gradient of logL with respect to
        theta.
        """
        raise NotImplementedError

    def information(self, theta, **extra):
        """
        Score function of model = - Hessian of logL with respect to
        theta.
        """
        raise NotImplementedError

    def newton(self, theta, **extra):
        def f(theta):
            return -self.logL(theta)
        self.results = scipy.optimize.fmin(f, theta)
        
class LikelihoodModelResults:

    def __init__(self, beta, normalized_cov_beta=None, scale=1.):
        self.beta = beta
        self.normalized_cov_beta = normalized_cov_beta
        self.scale = 1.
        self.sd = N.sqrt(self.scale)

    def t(self, column=None):
        """
        Return the t-statistic for a given parameter estimate.

        Use Tcontrast for more complicated t-statistics.

        """

        if self.normalized_cov_beta is None:
            raise ValueError, 'need covariance of parameters for computing T statistics'

        if column is None:
            _t = N.zeros(_beta.shape, N.float64)
            for i in range(self.beta.shape[0]):
                _t[i] = _beta[i] * recipr((self.sd * self.sqrt(self.normalized_cov_beta[i,i])))
        else:
            i = column
            _t = _beta[i] * recipr((self.sd * self.sqrt(self.normalized_cov_beta[i,i])))
        return _t

    def cov_beta(self, matrix=None, column=None, scale=None, other=None):
        """
        Returns the variance/covariance matrix of a linear contrast
        of the estimates of beta, multiplied by scale which
        will usually be an estimate of sigma^2.

        The covariance of
        interest is either specified as a (set of) column(s) or a matrix.
        """
        if self.normalized_cov_beta is None:
            raise ValueError, 'need covariance of parameters for computing (unnormalized) covariances'

        if scale is None:
            scale = self.scale

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
        
        if self.normalized_cov_beta is None:
            raise ValueError, 'need covariance of parameters for computing T statistics'

        results = ContrastResults()

        results.effect = N.dot(matrix, self.beta)
        if sd:
            results.sd = N.sqrt(self.cov_beta(matrix=matrix)) * self.sd
        if t:
            results.t = results.effect * recipr(results.sd)
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
        
        if self.normalized_cov_beta is None:
            raise ValueError, 'need covariance of parameters for computing F statistics'

        results = ContrastResults()
        cbeta = N.dot(matrix, self.beta)

        q = matrix.shape[0]
        if invcov is None:
            invcov = inv(self.cov_beta(matrix=matrix, scale=1.0))
        results.F = N.add.reduce(N.dot(invcov, cbeta) * cbeta, 0) * recipr((q * self.scale))
        return results

class ContrastResults:
    """
    Results from looking at a particular contrast of coefficients in
    a parametric model. The class does nothing, it is a container
    for the results from T and F contrasts.
    """
    pass

