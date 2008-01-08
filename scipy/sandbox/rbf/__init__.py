#
# rbf - Radial Basis Functions
#

from info import __doc__

from rbf import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from scipy.testing.pkgtester import Tester
test = Tester().test
