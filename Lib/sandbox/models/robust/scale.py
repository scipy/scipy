import numpy as N
import scipy
import scipy.stats

def MAD(a, c=0.6745):
    """
    Median Absolute Deviation along first axis of an array:

    median(abs(a - median(a))) / c

    """

    a = N.asarray(a, N.Float)
    d = N.multiply.outer(scipy.median(a), N.ones(a.shape[1:]))
    return scipy.median(N.fabs(a - d) / c)

class Huber:
    """
    Huber's proposal 2 for estimating scale.

    R Venables, B Ripley. \'Modern Applied Statistics in S\'
    Springer, New York, 2002.
    """

    c = 1.5
    tol = 1.0e-06

    tmp = 2 * scipy.stats.norm.cdf(c) - 1
    gamma = tmp + c**2 * (1 - tmp) - 2 * c * scipy.stats.norm.pdf(c)
    del(tmp)
    
    niter = 10

    def __call__(self, a, mu=None, scale=None):
        """
        Compute Huber\'s proposal 2 estimate of scale, using an optional
        initial value of scale and an optional estimate of mu. If mu
        is supplied, it is not reestimated.
        """

        self.a = N.asarray(a, N.Float)
        if mu is None:
            self.n = self.a.shape[0] - 1
            self.mu = N.multiply.outer(scipy.median(self.a), N.ones(self.a.shape[1:]))
            self.est_mu = True
        else:
            self.n = self.a.shape[0]
            self.mu = mu
            self.est_mu = False

        if scale is None:
            self.scale = MAD(self.a)**2
        else:
            self.scale = scale

        for donothing in self:
            pass

        self.s = N.sqrt(self.scale)
        return self.s

    def __iter__(self):
        self.iter = 0
        return self

    def next(self):
        a = self.a
        subset = self.subset(a)
        if self.est_mu:
            mu = (subset * a).sum() / a.shape[0]
        else:
            mu = self.mu
            
        scale = N.sum(subset * (a - mu)**2, axis=0) / (self.n * Huber.gamma - N.sum(1. - subset, axis=0) * Huber.c**2)

        self.iter += 1

        if (N.fabs(N.sqrt(scale) - N.sqrt(self.scale)) <= N.sqrt(self.scale) * Huber.tol and
            N.fabs(mu - self.mu) <= N.sqrt(self.scale) * Huber.tol):
            self.scale = scale
            self.mu = mu
            raise StopIteration
        else:
            self.scale = scale
            self.mu = mu

        if self.iter >= self.niter:
            raise StopIteration

    def subset(self, a):
        tmp = (a - self.mu) / N.sqrt(self.scale)
        return N.greater(tmp, -Huber.c) * N.less(tmp, Huber.c)

huber = Huber()    
