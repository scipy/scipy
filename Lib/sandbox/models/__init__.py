#
# models - Statistical Models
#

from scipy.sandbox.models.info import __doc__

import scipy.sandbox.models.model
import scipy.sandbox.models.formula
import scipy.sandbox.models.regression
import scipy.sandbox.models.robust
import scipy.sandbox.models.family
from scipy.sandbox.models.glm import model as glm
from scipy.sandbox.models.rlm import model as rlm

__all__ = filter(lambda s:not s.startswith('_'),dir())

from numpy.testing import NumpyTest
test = NumpyTest().test
