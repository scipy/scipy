#
# spatial - spatial data structures and algorithms
#

from info import __doc__

from kdtree import *


__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test

