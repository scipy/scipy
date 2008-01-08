#
# models - Statistical Models
#

__docformat__ = 'restructuredtext'

from scipy.stats.models.info import __doc__

import scipy.stats.models.model
import scipy.stats.models.formula
import scipy.stats.models.regression
import scipy.stats.models.robust
import scipy.stats.models.family
from scipy.stats.models.glm import Model as glm
from scipy.stats.models.rlm import Model as rlm

__all__ = filter(lambda s:not s.startswith('_'),dir())

from scipy.testing.pkgtester import Tester
test = Tester().test
