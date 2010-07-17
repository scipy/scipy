
# Modules contributed by BasSw (wegwerp@gmail.com)
from codata import *
from constants import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
