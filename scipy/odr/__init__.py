#
# odr - Orthogonal Distance Regression
#

from info import __doc__

__version__ = '0.7'
__author__ = 'Robert Kern <robert.kern@gmail.com>'
__date__ = '2006-09-21'

from odrpack import *
from models import *

__all__ = filter(lambda s: not s.startswith('_'), dir())

from numpy.testing import Tester
test = Tester().test
