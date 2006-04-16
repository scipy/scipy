"Linear Solvers"

from info import __doc__

from linsolve import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import ScipyTest
test = ScipyTest().test
