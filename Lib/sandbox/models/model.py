import scipy.optimize

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
        
