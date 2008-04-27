#
# cluster - Vector Quantization / Kmeans
#

from info import __doc__

__all__ = ['vq', 'hierarchy']

import vq, hierarchy
from scipy.testing.pkgtester import Tester
test = Tester().test
