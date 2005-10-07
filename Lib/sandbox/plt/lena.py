"""
lena
"""
__all__ = ['lena']

import os
import cPickle as pickle
from scipy_base import array

def lena():
    """ Return Lena picture as an array."""
    fname = os.path.join(os.path.dirname(__file__),'lena.dat')
    f = open(fname,'rb')
    arr = array(pickle.load(f))
    f.close()
    return arr

