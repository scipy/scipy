"""LOBPCG eigensolver"""

from info import __doc__
import lobpcg
if __doc__ and lobpcg.__doc__:
    __doc__ = '\n\n'.join( (lobpcg.__doc__, __doc__) )
del lobpcg

from lobpcg import *

from numpy.testing import Tester
test = Tester().test
