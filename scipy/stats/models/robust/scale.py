import numpy as N
from scipy.stats import norm, median

def unsqueeze(data, axis, oldshape):
    """
    unsqueeze a collapsed array

    >>> from numpy import mean
    >>> from numpy.random import standard_normal
    >>> x = standard_normal((3,4,5))
    >>> m = mean(x, axis=1)
    >>> m.shape
    (3, 5)
    >>> unsqueeze(m, 1, x.shape)
    >>> m.shape
    (3, 1, 5)
    >>>                       
    """
    
    newshape = list(oldshape)
    newshape[axis] = 1
    data.shape = newshape
    

def MAD(a, c=0.6745, axis=0):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    """

    a = N.asarray(a, N.float64)
    d = median(a, axis=axis)
    unsqueeze(d, axis, a.shape)

    return median(N.fabs(a - d) / c, axis=axis)

class Huber:
    """
    Huber's proposal 2 for estimating scale.

    R Venables, B Ripley. \'Modern Applied Statistics in S\'
    Springer, New York, 2002.
    """

    c = 1.5
    tol = 1.0e-06

    tmp = 2 * norm.cdf(c) - 1
    gamma = tmp + c**2 * (1 - tmp) - 2 * c * norm.pdf(c)
    del(tmp)
    
    niter = 30

    def __call__(self, a, mu=None, scale=None, axis=0):
        """
        Compute Huber\'s proposal 2 estimate of scale, using an optional
        initial value of scale and an optional estimate of mu. If mu
        is supplied, it is not reestimated.
        """

        self.axis = axis
        self.a = N.asarray(a, N.float64)
        if mu is None:
            self.n = self.a.shape[0] - 1
            self.mu = median(self.a, axis=axis)
            self.est_mu = True
        else:
            self.n = self.a.shape[0]
            self.mu = mu
            self.est_mu = False

        if scale is None:
            self.scale = MAD(self.a, axis=self.axis)**2
        else:
            self.scale = scale

        unsqueeze(self.scale, self.axis, self.a.shape)
        unsqueeze(self.mu, self.axis, self.a.shape)

        for donothing in self:
            pass

        self.s = N.squeeze(N.sqrt(self.scale))
        del(self.scale); del(self.mu); del(self.a)
        return self.s

    def __iter__(self):
        self.iter = 0
        return self

    def next(self):
        a = self.a
        subset = self.subset(a)
        if self.est_mu:
            mu = N.sum(subset * a + (1 - Huber.c) * subset, axis=self.axis) / a.shape[self.axis]
        else:
            mu = self.mu
        unsqueeze(mu, self.axis, self.a.shape)
            
        scale = N.sum(subset * (a - mu)**2, axis=self.axis) / (self.n * Huber.gamma - N.sum(1. - subset, axis=self.axis) * Huber.c**2)

        self.iter += 1

        if N.alltrue(N.less_equal(N.fabs(N.sqrt(scale) - N.sqrt(self.scale)), N.sqrt(self.scale) * Huber.tol)) and N.alltrue(N.less_equal(N.fabs(mu - self.mu), N.sqrt(self.scale) * Huber.tol)):
            self.scale = scale
            self.mu = mu
            raise StopIteration
        else:
            self.scale = scale
            self.mu = mu

        unsqueeze(self.scale, self.axis, self.a.shape)

        if self.iter >= self.niter:
            raise StopIteration

    def subset(self, a):
        tmp = (a - self.mu) / N.sqrt(self.scale)
        return N.greater(tmp, -Huber.c) * N.less(tmp, Huber.c)

huber = Huber()    
