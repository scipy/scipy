#
# stats - Statistical Functions
#

from info import __doc__

from stats import *
from distributions import *
from rv import *
from morestats import *
from kde import gaussian_kde

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
