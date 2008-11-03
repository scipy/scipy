#
# spatial - Distances
#

from info import __doc__
from kdtree import *
from ckdtree import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
__all__ += ['distance']

import distance
from numpy.testing import Tester
test = Tester().test
