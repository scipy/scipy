"Sparse eigenvalue solvers"

from info import __doc__

from arpack import *
from lobpcg import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
