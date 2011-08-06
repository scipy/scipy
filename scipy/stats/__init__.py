#
# stats - Statistical Functions
#

from info import __doc__

from stats import *
from distributions import *
from rv import *
from morestats import *
from kde import gaussian_kde
import mstats

#remove vonmises_cython from __all__, I don't know why it is included
__all__ = filter(lambda s:not (s.startswith('_') or s.endswith('cython')),dir())

from numpy.testing import Tester
test = Tester().test
