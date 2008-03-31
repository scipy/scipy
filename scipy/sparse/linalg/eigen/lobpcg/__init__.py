"""LOBPCG eigensolver"""

from info import __doc__
import lobpcg
__doc__ = '\n\n'.join( (lobpcg.__doc__, __doc__) )
del lobpcg

from lobpcg import *

from scipy.testing.pkgtester import Tester
test = Tester().test
