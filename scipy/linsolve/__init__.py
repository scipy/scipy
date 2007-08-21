"Linear Solvers"

from info import __doc__

import umfpack
__doc__ = '\n\n'.join( (__doc__,  umfpack.__doc__) )
del umfpack

from linsolve import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import NumpyTest
test = NumpyTest().test
