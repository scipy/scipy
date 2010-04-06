#
# signal - Signal Processing Tools
#

from info import __doc__

import sigtools
from waveforms import *
from bsplines import *
from filter_design import *
from ltisys import *
from windows import *
from signaltools import *
from wavelets import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
