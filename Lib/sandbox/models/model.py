import numpy as N
from numpy.linalg import inv
#from scipy import optimize

from scipy.sandbox.models.contrast import ContrastResults
from scipy.sandbox.models.utils import recipr

class Model:
    """
    A (predictive) statistical model. The class Model itself does nothing
    but lays out the methods expected of any subclass.
    """

    def __init__(self):
        pass

    def initialize(self):
        """
        Initialize (possibly re-initialize) a Model instance. For
        instance, the design matrix of a linear model may change
        and some things must be recomputed.
        """
        raise NotImplementedError

    def fit(self): 
        """
        Fit a model to data.
        """
        raise NotImplementedError

    def predict(self, design=None):
        """
        After a model has been fit, results are (assumed to be) stored
        in self.results, which itself should have a predict method.
        """
        self.results.predict(design)

    def view(self):
        """
        View results of a model.
        """
        raise NotImplementedError

class likelihood_model(Model):

    def logL(self, theta):
        """
        Log-likelihood of model.
        """
        raise NotImplementedError

    def score(self, theta):
        """
        Score function of model = gradient of logL with respect to
        theta.
        """
        raise NotImplementedError

    def information(self, theta):
        """
        Score function of model = - Hessian of logL with respect to
        theta.
        """
        raise NotImplementedError

    def newton(self, theta):
        raise NotImplementedError
#         def f(theta):
#             return -self.logL(theta)
#         self.results = optimize.fmin(f, theta)
        
class LikelihoodModelResults:

    def __init__(self, beta, normalized_cov_beta=None, scale=1.):
        self.beta = beta
        self.normalized_cov_beta = normalized_cov_beta
        self.scale = 1.

    def t(self, column=None):
        """
        Return the t-statistic for a given parameter estimate.

        Use Tcontrast for more complicated t-statistics.

        """

        if self.normalized_cov_beta is None:
            raise ValueError, 'need covariance of parameters for computing T statistics'

        if column is None:
            column = range(self.beta.shape[0])

        column = N.asarray(column)
        _beta = self.beta[column]
        _cov = self.cov_beta(column=column)
        if _cov.ndim == 2:
            _cov = N.diag(_cov)
        _t = _beta * recipr(N.sqrt(_cov))
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

        _t = _sd = None

        _effect = N.dot(matrix, self.beta)
        if sd:
            _sd = N.sqrt(self.cov_beta(matrix=matrix))
        if t:
            _t = _effect * recipr(_sd)
        return ContrastResults(effect=_effect, t=_t, sd=_sd, df_denom=self.df_resid)

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

        cbeta = N.dot(matrix, self.beta)

        q = matrix.shape[0]
        if invcov is None:
            invcov = inv(self.cov_beta(matrix=matrix, scale=1.0))
        F = N.add.reduce(N.dot(invcov, cbeta) * cbeta, 0) * recipr((q * self.scale))
        return ContrastResults(F=F, df_denom=self.df_resid, df_num=invcov.shape[0])


