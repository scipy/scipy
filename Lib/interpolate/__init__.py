#
# interpolate - Interpolation Tools
#

from info import __doc__

from interpolate import *
from fitpack import *

# New interface to fitpack library:
from fitpack2 import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from scipy.testing import ScipyTest 
test = ScipyTest().test
