""" Basic utilities for the GA package.
"""

from numpy import mean, std

from prng import prng


class GAError(Exception):
    """ Error from the GA code.
    """

def nop(x):
    """ Basic 'no-op' stub useful for interfaces which require a function.
    """
    return x

def flip_coin(p):
    """ Return True with probability p.
    """
    return (prng.random() < p)

class empty_class:
    """ Dummy class for cloning objects.
    """
    pass

def shallow_clone(item):
    """ Make a simple clone of an object.

    The attributes are not copied, just referenced.
    """
    new = empty_class()
    new.__class__ = item.__class__
    new.__dict__.update(item.__dict__)
    return new

def remove_NaN(z):
    """ Return an array with only finite (non-NaN, non-inf) values.
    """
    from numpy import isfinite
    return z[isfinite(z)]

def my_std(s):
    """ Standard deviation robust to NaNs and infs.
    """
    a = remove_NaN(s)
    if len(a) > 1:
        return std(a)
    else:
        return 0.

def my_mean(s):
    """ Mean robust to NaNs and infs.
    """
    a = remove_NaN(s)
    if len(a) > 0:
        return mean(a)
    else:
        return 0.
