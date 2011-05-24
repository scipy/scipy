from info import __doc__
from maxentropy import *
from scipy.misc import logsumexp

from numpy.testing import Tester
test = Tester().test

import warnings
warnings.warn("""
The scipy.maxentropy module is deprecated in scipy 0.10, and scheduled to be
removed in 0.11.

If you are using some of the functionality in this module and are of the
opinion that it should be kept or moved somewhere - or you are even interested
to maintain/improve this whole module - please ask on the scipy-dev mailing
list.

The logsumexp function has already been moved to scipy.misc.""",
DeprecationWarning)
