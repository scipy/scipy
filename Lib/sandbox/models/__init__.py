#
# models - Statistical Models
#

from info import __doc__

import model
import formula
import regression
import robust
import family
from glm import model as glm
from rlm import model as rlm

from numpy.testing import NumpyTest
test = NumpyTest().test
