""" This module implements the Always Stable Predictor from Ref. [1] that can
be used for predicting values in a (time) series, such as in [2-3].

Given a (time) series of a data source, this method can be used to predict the
value of data source at the next time step. This is implemented in the method
`always_stable_predictor`.

Usually this is implemented with a sliding window approach that removes the
oldest element in a list and adds a new element to the list. This is currently
implemented in the `Cache` class.

[1] J. Kolafa, J. Comput. Chem. 25, 335 (2004).
[2] T. D. Kuehne, M. Krack, F. Mohamed, and M. Parrinello,
    Phys. Rev. Lett. 98, 066401 (2007).
[3] T. Spura, H. Elgabarty, and T. D. Kuehne, Phys. Chem. Chem. Phys. (2015).
"""

import numpy as np
from scipy import misc

__all__ = ["always_stable_predictor", "Cache"]


def always_stable_predictor(data):
    """
    Run the Always Stable Predictor algorithm.

    Parameters
    ----------
    data: array_like
        The observable data that should be predicted to the next time step.

    Returns
    -------
    The estimated value following the existing values in data.
    """
    ret = np.zeros_like(data[0])

    k = len(data)
    for m in range(1, k+1):
        n_m = k-m
        prefactor = (-1.0)**(m+1) * m * misc.comb(2*k, k-m) / misc.comb(2*k-2, k-1)
        ret += prefactor*data[n_m]
    return ret


class Cache(object):
    """ Class to help with predicting new observables.

    If the item parameter is not passed to the methods, a default (unnamed)
    cache is used.
    """
    def __init__(self, cache):
        """
        Initialize Cache object.

        Parameters
        ----------
        cache : integer
            The number of elements to be kept in the cache.
        """
        self.cache = cache
        self._cached = {}

    def add(self, value, item="default"):
        """
        Adds a new value to the cache.

        Parameters
        ----------
        item : hashable object
            Which cache should be accessed.
        value : an object
            The object to be added to the cache.

        Raises
        ------
        RuntimeError when the cache was initialized to be 0.
        """
        if item not in self._cached:
            self._cached[item] = []
        if self.cache > len(self._cached[item]):
            self._cached[item].append(value)
        elif self.cache > 0:
            del self._cached[item][0]
            self._cached[item].append(value)
        elif self.cache != 0:
            raise RuntimeError("Cache: NOTHING CACHED for %s!" % item)

    def predict(self, item="default"):
        """
        Predicts the values in the cache to the next time step.

        Parameters
        ----------
        item : hashable object
            Which cache should be accessed.

        Returns
        -------
        The predicted value.
        """
        if item not in self._cached:
            return None
        if not len(self._cached[item]):
            return None
        return always_stable_predictor(self._cached[item])

    def old(self, item="default"):
        """
        Return the most recently added item in the cache.

        Parameters
        ----------
        item : hashable object
            Which cache should be accessed.

        Returns
        -------
        The most recently added item in the cache.
        """
        return self._cached[item][-1]
