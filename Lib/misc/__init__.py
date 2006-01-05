
from info import __doc__

import limits
from common import *
from helpmod import *

try:
    from pilutil import *
except ImportError:
    pass

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import ScipyTest 
test = ScipyTest().test
