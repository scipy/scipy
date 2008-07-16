"Sparse Linear Algebra routines"

from info import __doc__

from isolve import *
from dsolve import *
from interface import *
from eigen import *


__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
