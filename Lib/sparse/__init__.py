"Rudimentary sparse matrix class"

from info import __doc__

from sparse import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import ScipyTest 
test = ScipyTest().test
