# montecarlo.py: Routines for Monte Carlo simulation

# Copyright: Ed Schofield, 2005-2006
# License: BSD-style (see LICENSE.txt at root of scipy tree)

__author__ = "Ed Schofield"

from __future__ import division
import random, math, bisect, string
import numpy
import scipy
#from scipy.montecarlo.intsampler import intsampler
from scipy.sandbox.montecarlo.intsampler import intsampler


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
        self.labels = numpy.array(mydict.keys(), object)
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
                return (sample, scipy.log(sampleprobs))
 

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

