
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

