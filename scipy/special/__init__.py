#
# special - Special Functions
#

from info import __doc__, __docformat__
#from special_version import special_version as __version__

from basic import *
import specfun
import orthogonal
from orthogonal import *
from spfun_stats import multigammaln
from lambertw import lambertw

__all__ = filter(lambda s:not s.startswith('_'),dir())

from numpy.dual import register_func
register_func('i0',i0)
del register_func

from numpy.testing import Tester
test = Tester().test
