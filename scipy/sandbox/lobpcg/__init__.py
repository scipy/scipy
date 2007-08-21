"LOBPCG"

from info import __doc__
import lobpcg
__doc__ = '\n\n'.join( (lobpcg.__doc__, __doc__) )
del lobpcg

from lobpcg import *

from numpy.testing import NumpyTest
test = NumpyTest().test
