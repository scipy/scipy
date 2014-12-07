import sys
import warnings
import numbers

import numpy as np


class DeprecatedImport(object):
    """
    Deprecated import, with redirection + warning.

    Examples
    --------
    Suppose you previously had in some module::

        from foo import spam

    If this has to be deprecated, do::

        spam = DeprecatedImport("foo.spam", "baz")

    to redirect users to use "baz" module instead.

    """

    def __init__(self, old_module_name, new_module_name):
        self._old_name = old_module_name
        self._new_name = new_module_name
        __import__(self._new_name)
        self._mod = sys.modules[self._new_name]

    def __dir__(self):
        return dir(self._mod)

    def __getattr__(self, name):
        warnings.warn("Module %s is deprecated, use %s instead"
                      % (self._old_name, self._new_name),
                      DeprecationWarning)
        return getattr(self._mod, name)


# copy-pasted from scikit-learn utils/validation.py
def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance

    If seed is None (or np.random), return the RandomState singleton used 
    by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)
