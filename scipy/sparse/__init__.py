"Rudimentary sparse matrix class"

from info import __doc__

from sparse import *
from construct import *
from spfuncs import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import NumpyTest
test = NumpyTest().test
