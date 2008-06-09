#
# cluster - Vector Quantization / K-means / Hierarchical Clustering
#

from info import __doc__

__all__ = ['vq', 'hierarchy', 'distance']

import vq, hierarchy
from scipy.testing.pkgtester import Tester
test = Tester().test
