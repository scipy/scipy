"Multigrid Solvers"

from info import __doc__

from multilevel import multilevel_solver
from rs import ruge_stuben_solver
from sa import smoothed_aggregation_solver
from gallery import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from scipy.testing.pkgtester import Tester
test = Tester().test
