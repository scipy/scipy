
from info import __doc__

__all__ = ['limits']

import limits
from common import *
from helpmod import *

try:
    from pilutil import *
    __all__ += pilutil.__all__
except ImportError:
    pass

__all__ += common.__all__
__all__ += helpmod.__all__

from numpy.testing import ScipyTest
test = ScipyTest().test
