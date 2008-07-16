#
# spline - Spline Tools
#

from info import __doc__

from fitpack import *

# New interface to fitpack library:
from spline import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
