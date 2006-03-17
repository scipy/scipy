# montecarlo.py: Routines for Monte Carlo simulation

# Copyright: Ed Schofield, 2005-2006
# License: BSD-style (see LICENSE.txt at root of scipy tree)

__author__ = "Ed Schofield"

from __future__ import division
import random, math, bisect, string
import numpy
import scipy
#from scipy.montecarlo.intsampler import intsampler
from scipy.sandbox.montecarlo._intsampler import _intsampler


class intsampler(object):
    """A class that samples objects from a given discrete distribution.
    The distribution is defined on an integer-valued sample space and
    specified with a PMF as a list or array like this:

    >>> table = [10, 15, 20]

    representing this pmf:
        x       0       1       2
        p(x)    10/45   15/45   20/45

    The output will be something like:
    >>> sampler = intsampler(table)

    >> sampler.sample(10)
    array([c, b, b, b, b, b, c, b, b, b], dtype=object)

    The class is a thin wrapper around the '_intsampler' class, which
    uses the 5 table compressed lookup sampler of Marsaglia (2004).  It
    is very fast, and requires near-constant time in the size of the
    sample space.  It provides a friendlier interface than _intsampler,
    allowing use of an integer-valued PMF (unnormalized) and able to
    return pmf values of the sample points in addition to the sample
    itself.
    """
    def __init__(self, mytable):
        self.probs = numpy.array(mytable, float)
        s = self.probs.sum()
        if s > 0:
            self.probs /= s
        else:
            raise ValueError, "sum of table frequencies must be > 0"

        self.sampler =  _intsampler(self.probs)

    def sample(self, size, return_probs=0):
        """Generates a sample of the given size from the specified
        discrete distribution, optionally returning the probabilities
        under the distribution.

        The optional argument return_probs, 0 by default, has the
        following meaning:
            0: don't return pmf values at each sample point
            1: return pmf values at each sample point
            2: return log pmf values at each sample point
        """
        sample = self.sampler.sample(size)
        if return_probs == 0:
            return sample
        elif return_probs > 0:
            # More fancy indexing:
            sampleprobs = self.probs[sample]
            if return_probs == 1:
                return (sample, sampleprobs)
            elif return_probs == 2:
                return (sample, numpy.log(sampleprobs))



class dictsampler(object):
    """A class that samples objects from a given discrete distribution.
    The distribution is specified as a dictionary representing a PMF
    like this:
    >>> table = {}
    >>> table['a'] = 10
    >>> table['b'] = 15
    >>> table['c'] = 20

    The output will be something like:
    >>> sampler = dictsampler(table)

    >> sampler.sample(10)
    array([c, b, b, b, b, b, c, b, b, b], dtype=object)

    giving a sample from this distribution:
    x       'a'       'b'       'c'
    p(x)   10/180   150/180   20/180

    The function uses the constant-time 'intsampler' class, and should be
    very fast.
    """
    def __init__(self, mydict):
        # We can't use this:
        #   self.labels = numpy.array(mydict.keys(), object)
        # since numpy's construction of object arrays is dodgy.  Instead,
        # create an empty object array and fill it:
        self.labels = numpy.empty(len(mydict), dtype=object)
        for i, label in enumerate(mydict):
            self.labels[i] = label
        self.probs = numpy.array(mydict.values(), float)
        s = self.probs.sum()
        if s > 0:
            self.probs /= s
        else:
            raise ValueError, "sum of table frequencies must be > 0"

        self.sampler =  intsampler(self.probs)

    def sample(self, size, return_probs=0):
        """Generates a sample of the given size from the specified
        discrete distribution, optionally returning the probabilities
        under the distribution.

        The optional argument return_probs, 0 by default, has the
        following meaning:
            0: don't return pmf values at each sample point
            1: return pmf values at each sample point
            2: return log pmf values at each sample point
        """
        sampleindices = self.sampler.sample(size)
        # Fancy indexing with the object array of labels
        sample = self.labels[sampleindices]
        if return_probs == 0:
            return sample
        elif return_probs > 0:
            # More fancy indexing:
            sampleprobs = self.probs[sampleindices]
            if return_probs == 1:
                return (sample, sampleprobs)
            elif return_probs == 2:
                return (sample, numpy.log(sampleprobs))


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
