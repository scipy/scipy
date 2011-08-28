
from numpy import vectorize
from numpy.random import random_sample

__all__ = ['randwppf', 'randwcdf']

# XXX: Are these needed anymore?

#####################################
# General purpose continuous
######################################

def randwppf(ppf, args=(), size=None):
    """
    returns an array of randomly distributed integers of a distribution
    whose percent point function (inverse of the CDF or quantile function)
    is given.

    args is a tuple of extra arguments to the ppf function (i.e. shape,
    location, scale), and size is the size of the output.  Note the ppf
    function must accept an array of q values to compute over.

    """
    U = random_sample(size=size)
    return apply(ppf, (U,)+args)

def randwcdf(cdf, mean=1.0, args=(), size=None):
    """returns an array of randomly distributed integers of a distribution
    whose cumulative distribution function (CDF) is given.

    mean is the mean of the distribution (helps the solver).
    args is a tuple of extra arguments to the cdf function (i.e. shape,
    location, scale), and size is the size of the output.  Note the
    cdf function needs to accept a single value to compute over.
    """
    import scipy.optimize as optimize
    def _ppfopt(x, q, *nargs):
        newargs = (x,)+nargs
        return cdf(*newargs) - q

    def _ppf(q, *nargs):
        return optimize.fsolve(_ppfopt, mean, args=(q,)+nargs)

    _vppf = vectorize(_ppf)
    U = random_sample(size=size)
    return apply(_vppf,(U,)+args)
