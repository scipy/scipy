"Multigrid Solvers"

from info import __doc__

from multilevel import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import NumpyTest
test = NumpyTest().test
