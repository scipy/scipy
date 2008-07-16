
from info import __doc__

__all__ = ['who', 'source', 'info']

from common import *
from numpy import who, source, info as _info

import sys
def info(object=None,maxwidth=76,output=sys.stdout,toplevel='scipy'):
    return _info(object, maxwidth, output, toplevel)
info.__doc__ = _info.__doc__
del sys

try:
    from pilutil import *
    __all__ += pilutil.__all__
except ImportError:
    pass

__all__ += common.__all__

from numpy.testing import Tester
test = Tester().test
