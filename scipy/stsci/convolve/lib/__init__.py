import sys
from Convolve import *
import iraf_frame

__version__ = '2.0'

def test(level=1, verbosity=1):
    from numpy.testing import Tester
    return Tester().test(level,verbosity)
