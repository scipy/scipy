
from info import __doc__

import limits
from common import *
from helpmod import *
from pilutil import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from scipy.testing import ScipyTest 
test = ScipyTest().test
